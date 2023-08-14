import scipy.sparse as sparse
import random
import scipy.io as sio
import numpy as np
import multiprocessing
import time
import copy
import sys
import os

from sim_SortTime import spada_sortTime, feasta_sortTime, sparsecore_sortTime


wind_row = 8
wind_col = 8
processor_num = 128
hw_mode = "SparseCore"
assert hw_mode in ("FEASTA_1*8*8", "FEASTA_2*8*4", "SparseCore")
# if hw_mode == "SparseCore", feasta_time means sparsecore time
mtx_data_dir = "/home/eva_share/zhongkai/graph_set/"


def get_sparse_matrix_AB_csr(file_dir, f):
    sp_mat_A = sio.mmread(file_dir)
    sp_mat_B = copy.deepcopy(sp_mat_A)
    
    cols_A, rows_A = sp_mat_A.shape

    if cols_A != rows_A:
        sp_mat_B = sp_mat_B.transpose()
    
    sp_mat_A = sp_mat_A.tocsr()
    f.write("A nonzeros: " + str(sp_mat_A.nnz) + "\n")
    sp_mat_B = sp_mat_B.tocsr()
    f.write("B nonzeros: " + str(sp_mat_B.nnz) + "\n")

    assert sp_mat_A.has_sorted_indices
    assert sp_mat_B.has_sorted_indices

    cols_B, rows_B = sp_mat_B.shape
    f.write("A shape: " + str(sp_mat_A.shape) + "\n")
    f.write("B shape: " + str(sp_mat_B.shape) + "\n")
    # print(sp_mat_A.indptr.shape)
    # print(sp_mat_B.indptr.shape)
    assert cols_A == rows_B
    assert rows_A == cols_B
    return sp_mat_A, sp_mat_B

def simulate_time(col_block, sp_mat_B):
    time_spada, non_sum =  spada_sortTime(col_block, sp_mat_B)
    if hw_mode != "SparseCore":
        time_feasta = feasta_sortTime(col_block, sp_mat_B)
    else:
        time_feasta = sparsecore_sortTime(col_block, sp_mat_B)
    return time_spada, time_feasta, non_sum

def calc_single_tile(args_list):
    sp_mat_A, sp_mat_B, col_tile_id = args_list[0], args_list[1], args_list[2]
    # print("col tile: ", str(col_tile_id))
    cols_A, rows_A = sp_mat_A.shape
    if hw_mode in ("FEASTA_1*8*8", "FEASTA_2*8*4"):
        max_row_tile_num = (rows_A + wind_row - 1) // wind_row
    elif hw_mode == "SparseCore": # A的列不切分
        max_row_tile_num = 1
    else:
        raise NotImplementedError
    start_col_id = wind_col * col_tile_id
    end_col_id = min(cols_A, wind_col * col_tile_id + wind_col)
    col_num_in_tile = end_col_id - start_col_id
    spata_col_time = np.zeros((col_num_in_tile, max_row_tile_num))
    feasta_col_time = np.zeros((col_num_in_tile, max_row_tile_num))
    non_element = np.zeros((col_num_in_tile, max_row_tile_num))
    if hw_mode in ("FEASTA_1*8*8", "FEASTA_2*8*4"):
        for col_id_in_tile in range(col_num_in_tile):
            this_col = sp_mat_A.indices[sp_mat_A.indptr[col_id_in_tile + start_col_id] : sp_mat_A.indptr[col_id_in_tile + start_col_id + 1]]
            row_tile_num = (len(this_col) + wind_row - 1) // wind_row
            for row_tile_id in range(row_tile_num):
                row_start = row_tile_id * wind_row
                row_end = min(row_tile_id * wind_row + wind_row, len(this_col))
                col_block = this_col[row_start : row_end]
                # print(col_block)
                spata_col_time[col_id_in_tile, row_tile_id], feasta_col_time[col_id_in_tile, row_tile_id], non_element[col_id_in_tile, row_tile_id] = simulate_time(col_block, sp_mat_B)
    elif hw_mode == "SparseCore":
        for col_id_in_tile in range(col_num_in_tile):
            this_col = sp_mat_A.indices[sp_mat_A.indptr[col_id_in_tile + start_col_id] : sp_mat_A.indptr[col_id_in_tile + start_col_id + 1]]
            spata_col_time[col_id_in_tile, 0], feasta_col_time[col_id_in_tile, 0], non_element[col_id_in_tile, 0] = simulate_time(this_col, sp_mat_B)
            spata_col_time[col_id_in_tile, 0] = 0
    else:
        raise NotImplementedError


    def col_time_reduction(col_time, reduction_ratio, total_cols_A):
        col_time_after_reduction = np.zeros((col_num_in_tile // reduction_ratio, max_row_tile_num)) # A的reduction_ratio行取max再累加
        for col_reduction_id in range(col_num_in_tile // reduction_ratio):
            col_reduction_start = col_reduction_id * reduction_ratio
            col_reduction_end = min(col_reduction_id * reduction_ratio + reduction_ratio, total_cols_A)
            col_time_after_reduction[col_reduction_id, :] = np.max(col_time[col_reduction_start : col_reduction_end, :], 0)
        return col_time_after_reduction


    if hw_mode == "FEASTA_2*8*4":
        feasta_col_time_reduction = col_time_reduction(feasta_col_time, 2, cols_A) # 第一维并行度为2，窗内A的每2行时间先取max再累加
    elif hw_mode == "FEASTA_1*8*8":
        feasta_col_time_reduction = feasta_col_time
    elif hw_mode == "SparseCore":
        feasta_col_time_reduction = col_time_reduction(feasta_col_time, 8, cols_A) # 第一维并行度为8，窗内A的每8行时间先取max再累加
        assert feasta_col_time_reduction.shape[1] == 1
    else:
        raise NotImplementedError

    spata_tile_time = np.max(spata_col_time, 0)
    feasta_tile_time = np.sum(feasta_col_time_reduction, 0)
    assert spata_tile_time.shape == (max_row_tile_num,)
    assert feasta_tile_time.shape == (max_row_tile_num,)
    spata_all_rows_tile_time = np.sum(spata_tile_time)
    feasta_all_rows_tile_time = np.sum(feasta_tile_time)
    non_element_sum = np.sum(non_element)
    # print("spata time: ", str(spata_all_rows_tile_time), "feasta time: ", str(feasta_all_rows_tile_time), "non_element: ", str(non_element_sum))
    
    return (col_tile_id, spata_all_rows_tile_time, feasta_all_rows_tile_time, non_element_sum)



def simulation_top(file_dir, output_dir):
    with open(output_dir, "a") as f:
        f.write("==============================================\n")
        f.write("%s: %s\n" % (file_dir, time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())))
        print("%s: %s" % (file_dir, time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())))
        print("loading mtx file: " + str(file_dir))
        sp_mat_A, sp_mat_B = get_sparse_matrix_AB_csr(file_dir, f)
    
        cols_A, rows_A = sp_mat_A.shape
        max_col_tile_num = (cols_A + wind_col - 1) // wind_col

        results = []
        multiprocess_time = time.time()

        print("preparing args, max_col_tile_num = ", str(max_col_tile_num))
        args_list = []
        for col_tile_id in range(max_col_tile_num):
            args = (sp_mat_A, sp_mat_B, col_tile_id)
            args_list.append(args)

        multiprocess_time = time.time()
        if True:
            ## 将参数传递给线程池，绑定执行方法，map方法返回的是一个结果列表，包含各个线程的执行结果
            print("multi processing")
            pool = multiprocessing.Pool(processes = processor_num)
            results = pool.map(calc_single_tile, args_list)
        else:
            print("single processing")
            for args in args_list:
                results.append(calc_single_tile(args))
    
        multiprocess_time = time.time() - multiprocess_time
        print("processing time: %.2f s" % multiprocess_time)

        print("post processing")
    
        spada_time = 0
        feasta_time = 0
        non_element = 0
        for col_tile_id in range(max_col_tile_num):
            assert col_tile_id == results[col_tile_id][0]
            spada_time += results[col_tile_id][1]
            feasta_time += results[col_tile_id][2]
            non_element += results[col_tile_id][3]


        f.write("spada  cycle: " + str(spada_time) + "\n")
        f.write("feasta cycle: " + str(feasta_time) + "\n")
        f.write("non_element: " + str(non_element) + "\n")
        print("spada  cycle: " + str(spada_time))
        print("feasta cycle: " + str(feasta_time))
        print("non_element: " + str(non_element))


if __name__ == "__main__":
    base_dir = "/home/eva_share/zhongkai/tf_set/"
    for name in ["bert-base-uncased-l1-query-0.9500.mtx", "opt-6.7b-attention-map-0.9813-7-22.mtx"]:
        simulation_top(os.path.join(base_dir, name), "result_tf_spmm.txt")

