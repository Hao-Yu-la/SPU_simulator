
MAXNUM = 1000000000
spmm_dim = 1024
cal_mode = "spmm"
assert cal_mode in ("spmspm", "spmm", "spmm_fix")
hw_mode = "FEASTA_1*8*8"
assert hw_mode in ("FEASTA_1*8*8", "FEASTA_2*8*4")

def spada_sortTime(a_1_8, b):
    # a_1_8: 1x8 list
    # b: matrix
    # return: time

    if cal_mode == "spmspm":
        b_line = [] # 与a中一行里的元素对应相乘的b中的行的列坐标
        non_sum = 0 # a中一行里的元素对应相乘的b中的行的非零元素个数
        for i in range(len(a_1_8)):
            # print("a_1_8[i]: ", a_1_8[i])
            # print(b.indices.shape)
            # print(b.indptr.shape)
            b_line.append(list(b.indices[b.indptr[a_1_8[i]] : b.indptr[a_1_8[i] + 1]]))
            non_sum += len(b_line[i])

            # print("b_line: ", b_line)
        time_line = spada_sortTime_1D(b_line) # 计算c一行的时间（a一行中的元素乘b对应的行，然后merge）
    elif cal_mode == "spmm":
        b_col = b.shape[1]
        time_line = b_col
        non_sum = b_col * len(a_1_8)
    elif cal_mode == "spmm_fix":
        b_col = spmm_dim
        time_line = b_col
        non_sum = b_col * len(a_1_8)
    else:
        print("Error: cal_mode is wrong!")
        exit(1)

    return time_line, non_sum

def spada_sortTime_1D(b_line_in):
    # b_line: 2D list
    # return: time

    time = 0
    line_num = len(b_line_in)
    # print("b_line_in: ", b_line_in)
    b_line = [[] for _ in range(line_num)] 
    # print("b_line: ", b_line)

    # b中对应的行进行sort&merge的时间
    while b_line != [[]] * line_num or b_line_in != [[]] * line_num: 
        # b_line_in中的元素出队加入b_line
        # for i in range(line_num):
        #     if b_line_in[i] != []:
        #         b_line[i].append(b_line_in[i].pop(0))
        i = 0
        while i < line_num:
            if i == line_num - 1:
                if b_line_in[i] != []:
                    b_line[i].append(b_line_in[i].pop(0))
                i += 1
            else:
                if len(b_line_in[i]) > 1 and len(b_line_in[i+1]) > 0 and b_line_in[i][1] < b_line_in[i+1][0]:
                    b_line[i].append(b_line_in[i].pop(0))
                    b_line[i].append(b_line_in[i].pop(0))
                elif len(b_line_in[i]) > 0 and len(b_line_in[i+1]) > 1 and b_line_in[i][0] > b_line_in[i+1][1]:
                    b_line[i+1].append(b_line_in[i+1].pop(0))
                    b_line[i+1].append(b_line_in[i+1].pop(0))
                else:
                    if b_line_in[i] != []:
                        b_line[i].append(b_line_in[i].pop(0))
                    if b_line_in[i+1] != []:
                        b_line[i+1].append(b_line_in[i+1].pop(0))
                i += 2          
        # print("b_line_in: ", b_line_in)
        # print("b_line: ", b_line)

         # 找出每行第3个元素中的最小值
        min_element = MAXNUM
        min_index = -1
        for i in range(line_num):
            if b_line[i] == []:
                continue
            if len(b_line[i]) < 3:
                if b_line[i][len(b_line[i])-1] < min_element:
                    min_element = b_line[i][len(b_line[i])-1]
                    min_index = i
            else:
                if b_line[i][2] < min_element:
                    min_element = b_line[i][2]
                    min_index = i
            # print("min_element: ", min_element)

        # 每行前2个元素中小于min_element的元素出队
        # for i in range(line_num): 
        #     while b_line[i] != [] and b_line[i][0] < min_element:
        #         b_line[i].pop(0)
        # if b_line_in[min_index] == []: # 特殊情况 b_line_in[min_index]为空
        #     b_line[min_index].pop(0)
        #     time += 1
        for i in range(2):
            for j in range(line_num):
                if b_line[j] != [] and b_line[j][0] <= min_element:
                    b_line[j].pop(0)
        # print(min_element)
        # print("b_line: ", b_line)
        
        time += 1

    return time

def feasta_sortTime(a_1_8, b):
    # a_1_8: 1x8 list
    # b: matrix
    # return: time

    if cal_mode == "spmspm":
        b_line = [] # 与a中一行里的元素对应相乘的b中的行的列坐标
        for i in range(len(a_1_8)):
            b_line.append(list(b.indices[b.indptr[a_1_8[i]] : b.indptr[a_1_8[i] + 1]]))
        time_line = feasta_sortTime_1D(b_line) # 计算c一行的时间（a一行中的元素乘b对应的行，然后merge）
    elif cal_mode == "spmm":
        b_col = b.shape[1]
        if hw_mode == "FEASTA_2*8*4":
            time_line = b_col // 4 + 1
        else:
            time_line = b_col // 8 + 1
    elif cal_mode == "spmm_fix":
        b_col = spmm_dim
        if hw_mode == "FEASTA_2*8*4":
            time_line = b_col // 4 + 1
        else:
            time_line = b_col // 8 + 1
    else:
        print("Error: cal_mode is wrong!")
        exit(1)

    return time_line

def feasta_sortTime_1D(b_line):
    # b_line: 2D list
    # return: time

    time = 0
    line_num = len(b_line)

    # b中对应的行进行sort&merge的时间
    while b_line != [[]] * line_num : 

         # 找出每行第9/17个元素中的最小值
        if hw_mode == "FEASTA_2*8*4":
            look_num = 8
        else:
            look_num = 16
        min_element = MAXNUM
        min_index = -1
        for i in range(line_num):
            if b_line[i] == []:
                continue
            if len(b_line[i]) < (look_num + 1):
                if b_line[i][len(b_line[i])-1] < min_element:
                    min_element = b_line[i][len(b_line[i])-1]
                    min_index = i
            else:
                if b_line[i][look_num] < min_element:
                    min_element = b_line[i][look_num]
                    min_index = i

        # 输出小于min_element的元素，若多于64个则输出最小的64个
        if hw_mode == "FEASTA_2*8*4":
            pop_num = 32
        else:
            pop_num = 64
        for i in range(pop_num):
            min_element_0 = min_element
            min_index_0 = -1
            for j in range(line_num):
                if b_line[j] == []:
                    continue
                if b_line[j][0] <= min_element_0:
                    min_element_0 = b_line[j][0]
                    min_index_0 = j
            if min_index_0 != -1:
                b_line[min_index_0].pop(0)
            else:
                break
        # if b_line[min_index] == [min_element]:
        #     b_line[min_index].pop(0)
        #     time += 1

        time += 1

    return time

def sparsecore_sortTime(a_1_8, b):
    # a_1_8: 1x8 list
    # b: matrix
    # return: time

    if cal_mode == "spmspm":
        b_line = [] # 与a中一行里的元素对应相乘的b中的行的列坐标
        for i in range(len(a_1_8)):
            b_line.append(list(b.indices[b.indptr[a_1_8[i]] : b.indptr[a_1_8[i] + 1]]))
        if len(b_line) == 0:
            time_line = 0
        elif len(b_line) == 1:
            time_line = len(b_line[0]) // 8 + 1
        else:
            time_line = sparsecore_sortTime_1D(b_line) # 计算c一行的时间（a一行中的元素乘b对应的行，然后逐个merge，每次merge两行）
    elif cal_mode == "spmm":
        b_col = b.shape[1]
        time_line = b_col // 8 + 1 + (b_col // 4 + 1) * (len(a_1_8) - 1)
    elif cal_mode == "spmm_fix":
        b_col = spmm_dim
        time_line = b_col // 8 + 1 + (b_col // 4 + 1) * (len(a_1_8) - 1)
    else:
        print("Error: cal_mode is wrong!")
        exit(1)

    return time_line

def sparsecore_sortTime_1D(b_line):
    # b_line: 2D list
    # return: time

    time_sum = 0
    merge_line = b_line[0][:]
    line_num = len(b_line)
    for i in range(1,line_num):
        merge_line, time = sparsecore_merge(merge_line, b_line[i])
        time_sum += time

    return time_sum

def sparsecore_merge(merge_line, b_line_one):
    # merge_line: 1D list
    # b_line_one: 1D list
    # return: merge_line, time

    merge_line_new = []
    time = 0
    look_num = 8
    while merge_line != [] and b_line_one != []:
        if len(merge_line) < (look_num + 1) and len(b_line_one) < (look_num + 1):
            min_element = min(merge_line[-1], b_line_one[-1])
        elif len(merge_line) < (look_num + 1):
            min_element = min(merge_line[-1], b_line_one[look_num])
        elif len(b_line_one) < (look_num + 1):
            min_element = min(merge_line[look_num], b_line_one[-1])
        else:
            min_element = min(merge_line[look_num], b_line_one[look_num])
        
        # 输出小于min_element的元素，若多于8个则输出最小的8个
        for i in range(8):
            if merge_line != [] and b_line_one != []:
                if merge_line[0] <= min_element and b_line_one[0] <= min_element:
                    if merge_line[0] <= b_line_one[0]:
                        ele = merge_line.pop(0)
                        if merge_line_new == [] or ele != merge_line_new[-1]:
                            merge_line_new.append(ele)
                    else:
                        ele = b_line_one.pop(0)
                        if merge_line_new == [] or ele != merge_line_new[-1]:
                            merge_line_new.append(ele)
                elif merge_line[0] <= min_element:
                    ele = merge_line.pop(0)
                    if merge_line_new == [] or ele != merge_line_new[-1]:
                        merge_line_new.append(ele)
                elif b_line_one[0] <= min_element:
                    ele = b_line_one.pop(0)
                    if merge_line_new == [] or ele != merge_line_new[-1]:
                        merge_line_new.append(ele)
                else:
                    break
            elif merge_line != []:
                if merge_line[0] <= min_element:
                    ele = merge_line.pop(0)
                    if merge_line_new == [] or ele != merge_line_new[-1]:
                        merge_line_new.append(ele)
                else:
                    break
            elif b_line_one != []:
                if b_line_one[0] <= min_element:
                    ele = b_line_one.pop(0)
                    if merge_line_new == [] or ele != merge_line_new[-1]:
                        merge_line_new.append(ele)
                else:
                    break
            else:
                break
        time += 1
    
    if merge_line != []:
        time += len(merge_line)//8 + 1
        if merge_line_new == [] or merge_line_new[-1] != merge_line[0]:
            merge_line_new.extend(merge_line)
        else:
            merge_line_new.extend(merge_line[1:])
    elif b_line_one != []:
        time += len(b_line_one)//8 + 1
        if merge_line_new == [] or merge_line_new[-1] != b_line_one[0]:
            merge_line_new.extend(b_line_one)
        else:
            merge_line_new.extend(b_line_one[1:])
    
    return merge_line_new, time


        


if __name__ == '__main__':
    a = [0,1,2,3,4,5,6,7]
    b = [[2,5,9,11,12,14,15,16,17,18,19], [1,5,8,15,17,19,20,21,22,23,24], [3,6,7,10,11,12,14,16,18,19,21,22], [4,5,6,10,11,12,13,14,15,16,17,18], [1,2,3,5,9,10,15,17,18,19,20], [4,5,6,8,9,10,11,12,13,14,15,16], [7,8,9,10,11,15,19,20,21,22,25,26], [1,2,3,7,8,9,10,11,12,13,14,15]]
    print(spada_sortTime(a, b))
    print(a, b)
    print(feasta_sortTime(a, b))
    

