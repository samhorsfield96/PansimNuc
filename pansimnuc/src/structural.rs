// in this script, genes can move around, be duplicated and deleted

use crate::mutation::MutationMap;
use std::collections::HashMap;

// TODO need to think of way of rearranging each feature, and taking into account where insertions
// and translocations occur. Also need two different TE compartments, one which copies and inserts
// and one which moves around and inserts. 

// also need to think about how to determine whether a TE inserts into another gene, making it non-functional
// or whether it is upstream or downstream and can augment its function, having a multiplicative effect on its fitness contribution.