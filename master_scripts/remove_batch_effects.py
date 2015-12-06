# This script takes in a dataset and corresponding sample table and
# calls two big functions:
# 1) batch_correction.main(), which performs multiclass LDA based on
#    the batches present in the specified sample table column; and
# 2) check_batch_effects.main(), which takes in the stacked datasets
#    exported by the batch_correction script and visualizes the
#    change in 


import batch_correction
import check_batch_effects
