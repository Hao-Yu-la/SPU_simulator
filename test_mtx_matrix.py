import scipy.sparse as sparse
import random
import scipy.io as sio
import numpy as np
import multiprocessing
import time
import copy
from sim_SortTime import spada_sortTime, feasta_sortTime


wind_row = 8
wind_col = 8
processor_num = 128
mtx_data_dir = "/home/eva_share/zhongkai/graph_set/web-Google.mtx"

# np.random.seed(42)
# sparse_matrix = sparse.random(500, # colomn
#                 5000, # row
#                 density=0.05,
#                 format="csr")
# sio.mmwrite("sparse_matrix.mtx", sparse_matrix)
# print(sparse_matrix)


def get_sparse_matrix_AB_csr(file_dir):
    sp_mat_A = sio.mmread(file_dir)
    sp_mat_B = copy.deepcopy(sp_mat_A)
    
    cols_A, rows_A = sp_mat_A.shape

    if cols_A != rows_A:
        sp_mat_B = sp_mat_B.transpose()
    
    sp_mat_A = sp_mat_A.tocsr()
    print("A nonzeros: ", str(sp_mat_A.nnz))
    sp_mat_B = sp_mat_B.tocsr()
    print("B nonzeros: ", str(sp_mat_B.nnz))

    assert sp_mat_A.has_sorted_indices
    assert sp_mat_B.has_sorted_indices

    cols_B, rows_B = sp_mat_B.shape
    print("A shape: ", str(sp_mat_A.shape))
    print("B shape: ", str(sp_mat_B.shape))
    # print(sp_mat_A.indptr.shape)
    # print(sp_mat_B.indptr.shape)
    assert cols_A == rows_B
    assert rows_A == cols_B
    return sp_mat_A, sp_mat_B

def simulate_time(col_block, sp_mat_B):
    time_spada, non_sum =  spada_sortTime(col_block, sp_mat_B)
    time_feasta = feasta_sortTime(col_block, sp_mat_B)
    return time_spada, time_feasta, non_sum

def calc_single_tile(args_list):
    sp_mat_A, sp_mat_B, col_tile_id = args_list[0], args_list[1], args_list[2]
    print("col tile: ", str(col_tile_id))
    cols_A, rows_A = sp_mat_A.shape
    max_row_tile_num = (rows_A + wind_row - 1) // wind_row
    start_col_id = wind_col * col_tile_id
    end_col_id = min(cols_A, wind_col * col_tile_id + wind_col)
    col_num_in_tile = end_col_id - start_col_id
    spata_col_time = np.zeros((col_num_in_tile, max_row_tile_num))
    feasta_col_time = np.zeros((col_num_in_tile, max_row_tile_num))
    non_element = np.zeros((col_num_in_tile, max_row_tile_num))
    for col_id_in_tile in range(col_num_in_tile):
        this_col = sp_mat_A.indices[sp_mat_A.indptr[col_id_in_tile + start_col_id] : sp_mat_A.indptr[col_id_in_tile + start_col_id + 1]]
        row_tile_num = (len(this_col) + wind_row - 1) // wind_row
        for row_tile_id in range(row_tile_num):
            row_start = row_tile_id * wind_row
            row_end = min(row_tile_id * wind_row + wind_row, len(this_col))
            col_block = this_col[row_start : row_end]
            # print(col_block)
            spata_col_time[col_id_in_tile, row_tile_id], feasta_col_time[col_id_in_tile, row_tile_id], non_element[col_id_in_tile, row_tile_id] = simulate_time(col_block, sp_mat_B)
    
    spata_tile_time = np.max(spata_col_time, 0)
    feasta_tile_time = np.sum(feasta_col_time, 0)
    assert spata_tile_time.shape == (max_row_tile_num,)
    assert feasta_tile_time.shape == (max_row_tile_num,)
    spata_all_rows_tile_time = np.sum(spata_tile_time)
    feasta_all_rows_tile_time = np.sum(feasta_tile_time)
    non_element_sum = np.sum(non_element)
    # print("spata time: ", str(spata_all_rows_tile_time), "feasta time: ", str(feasta_all_rows_tile_time), "non_element: ", str(non_element_sum))
    
    return (col_tile_id, spata_all_rows_tile_time, feasta_all_rows_tile_time, non_element_sum)


if __name__ == "__main__":
    
    print("loading mtx")
    sp_mat_A, sp_mat_B = get_sparse_matrix_AB_csr(mtx_data_dir)
    
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

    print("spada  cycle: ", str(spada_time))
    print("feasta cycle: ", str(feasta_time))
    print("non_element: ", str(non_element))
