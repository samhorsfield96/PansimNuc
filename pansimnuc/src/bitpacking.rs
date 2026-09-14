// idea is to have two bits per base in a standard u64, and a separate u64 vector 
// holding positions of any Ns in a sequence
// main issue is dealing with indels, need to have a way of deleting element (could have positional deletion or addition?)

// seq holds the two bit characters of each base
// pos_n holds the one bit character of each N base if present
pub struct bitpacked {
    pos_base: u64
    pos_n: u32
}

impl bitpacked {
    // new function should take a slice of string and fill the bitpacked item
    fn new()

    // insert function should enable editing of a specific location within the bitpacked item
    fn insert()
}

// this is a full element, for which the length needs to be determined for account for the correct
// number of bitpacked elements
pub bitpacked_vec {
    seq: Vec<bitpacked>
}

impl bitpacked_vec {
    // new function should take a full string, create a series of bitpacked elements and fill them using the above functions
    fn new()

    // insert function should enable editing of a specific location within the bitpacked item
    fn insert()
}


