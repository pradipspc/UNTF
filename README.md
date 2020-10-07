We are providing codes for generating unit norm tight frames (UNTFs) with low coherence from biadjacency matrices asociated with bi-regular graphs via an existing embedding operation. Next, we show the sparse recovery performance of such UNTFs.  The biadjacency matrices asociated with bi-regular graphs have same number of ones in each column and row.
untf_mod.m generates UNTFs from binary matrices having constant row and column weight.
uob_binarymats.m generates binary matrices by evaluating polynomials over finite filelds. 
test_Bin_Gauss_with_OMP.m and test_UNTF_Gauss_with_BP.m evaluate sparse recovery performance of binary matrices via OMP and sparse recovery performance of UNTFs via BP, respectively.
