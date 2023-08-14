# SPU_simulator
|- README.md 

|- test_graph_matrix.py 

|- test_tf_matrix.py

|- test_cnn_matrix.py

|- sim_SortTime.py

- test_graph_matrix.py: test graph datasets
- test_tf_matrix.py: test transformer datasets
- test_cnn_matrix.py: test cnn datasets
- sim_SortTime.py: simulate the sorting merge time of the matrix multiplication

In test.py file, "hw_mode" is the mode of the hardware, "hw_mode" = FEASTA_1\*8\*8 means the parallel is 1\*8\*8, "hw_mode" = FEASTA_2\*8\*4 means the parallel is 2\*8\*4, "hw_mode" = SparseCore means the parallel is 8\*1\*8.\

In sim_SortTime.py, "hw_mode" is the mode of the hardware, "hw_mode" = FEASTA_1\*8\*8 means the parallel is 1\*8\*8, "hw_mode" = FEASTA_2\*8\*4 means the parallel is 2\*8\*4, "hw_mode" = SparseCore means the parallel is 8\*1\*8. But the "hw_mode" = SparseCore is not used in sim_SortTime.py, because when the parallel is 8\*1\*8, the test.py uses the sparsecore_sortTime function directly.
"cal_mode" is the mode of the calculation, "cal_mode" = spmspm means the calculation is the sparse matrix times sparse matrix, "cal_mode" = spmm means the calculation is the sparse matrix times dense matrix and the dense matrix column is decided by the dense matrix in the datasets, "cal_mode" = spmm_fix means the calculation is the sparse matrix times dense matrix and the dense matrix column is fixed by the parameter "spmm_dim" in the sim_SortTime.py file.

Change the "hw_mode" and "cal_mode" to test the different hardware and calculation.

Three test.py files are used to test the different datasets.