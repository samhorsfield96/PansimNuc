use bitvec::prelude::*;

// idea is to have two bits per base in a standard u64, and a separate u64 vector 
// holding positions of any Ns in a sequence
// main issue is dealing with indels, need to have a way of deleting element (could have positional deletion or addition?)

pub type DnaBits = BitVec<u64, Lsb0>;

// seq holds the two bit characters of each base
// pos_n holds the one bit character of each N base if present
#[derive(Clone, PartialEq, Eq, Debug)]
pub struct bitpacked {
    bits : DnaBits
}

impl bitpacked {
    // new function should take a slice of string and fill the bitpacked item
    pub fn new(seq: &str) -> Self {
        let mut bits = BitVec::<u64, Lsb0>::with_capacity(seq.len() * 3);

        for base in seq.bytes() {
            let value = match base {
                b'A' => 0, // 000
                b'C' => 1, // 001
                b'G' => 2, // 010
                b'T' => 3, // 011
                _ => 4,    // 100
            };

            bits.push(value & 0b001 != 0);
            bits.push(value & 0b010 != 0);
            bits.push(value & 0b100 != 0);
        }

        Self {
            bits: bits
        }
    }

    pub fn new_vec(seq: Vec<u8>) -> Self {
        let mut bits = BitVec::<u64, Lsb0>::with_capacity(seq.len() * 3);

        for value in seq.iter() {
            bits.push(value & 0b001 != 0);
            bits.push(value & 0b010 != 0);
            bits.push(value & 0b100 != 0);
        }

        Self {
            bits: bits
        }
    }

    pub fn decode_int(&self) -> Vec<u8> {
        let sequence_length = self.bits.len() / 3;

        let mut seq = vec![0u8; sequence_length];

        for site in 0..sequence_length {
            let offset = site * 3;

            let value =
                (self.bits[offset] as u8)
                | ((self.bits[offset + 1] as u8) << 1)
                | ((self.bits[offset + 2] as u8) << 2);

            seq[site] = value;
        }

        seq
    }

    // function to generate string
    pub fn decode_dna(&self) -> Vec<u8> {
        let sequence_length = self.bits.len() / 3;

        let mut seq = vec![0u8; sequence_length];

        for site in 0..sequence_length {
            let offset = site * 3;

            let value =
                (self.bits[offset] as u8)
                | ((self.bits[offset + 1] as u8) << 1)
                | ((self.bits[offset + 2] as u8) << 2);

            let base = match value {
                0 => b'A',
                1 => b'C',
                2 => b'G',
                3 => b'T',
                _ => b'N',
            };

            seq[site] = base;
        }

        seq
    }

    pub fn initialise() -> Self {
        let bits = BitVec::<u64, Lsb0>::with_capacity(1 * 3);

        Self {
            bits: bits
        }
    }

    pub fn len(&self) -> usize {
        let sequence_length = self.bits.len() / 3;
        sequence_length
    }

    pub fn index(&self, site: usize) -> u8 {

        let offset = site * 3;

        (self.bits[offset] as u8)
            | ((self.bits[offset + 1] as u8) << 1)
            | ((self.bits[offset + 2] as u8) << 2)
    }

    pub fn update(&mut self, site: usize, new_value: u8) {
        assert!(site < self.len(), "site out of bounds");
        assert!(new_value <= 3, "DNA value must be between 0 and 3, 4 is N");

        let offset = site * 3;

        self.bits.set(offset, new_value & 0b001 != 0);
        self.bits.set(offset + 1, new_value & 0b010 != 0);
        self.bits.set(offset + 2, new_value & 0b100 != 0);
    }

    // insert function should enable editing of a specific location within the bitpacked item
    pub fn insert(&mut self, site: usize, new_value: u8) {
        assert!(site <= self.len(), "site out of bounds");
        let offset = site * 3;

        self.bits.insert(offset, new_value & 0b001 != 0);
        self.bits.insert(offset + 1, new_value & 0b010 != 0);
        self.bits.insert(offset + 2, new_value & 0b100 != 0);
    }

    pub fn remove(&mut self, site: usize) {
        assert!(site < self.len(), "site out of bounds");

        let offset = site * 3;

        self.bits.remove(offset + 2);
        self.bits.remove(offset + 1);
        self.bits.remove(offset);
    }

    pub fn iter(&self) -> impl DoubleEndedIterator<Item = u8> + '_ {
        (0..self.len()).map(|site| self.index(site))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_dna_sequence() {
        let sequence = bitpacked::new("ACGTN");

        assert_eq!(sequence.len(), 5);
        assert_eq!(sequence.decode_int(), vec![0, 1, 2, 3, 4]);
        assert_eq!(sequence.decode_dna(), b"ACGTN");
    }

    #[test]
    fn indexes_bases() {
        let sequence = bitpacked::new("ACGTN");

        assert_eq!(sequence.index(0), 0);
        assert_eq!(sequence.index(1), 1);
        assert_eq!(sequence.index(2), 2);
        assert_eq!(sequence.index(3), 3);
        assert_eq!(sequence.index(4), 4);
    }

    #[test]
    fn updates_base() {
        let mut sequence = bitpacked::new("AAAA");

        sequence.update(0, 1); // C
        sequence.update(1, 2); // G
        sequence.update(2, 3); // T

        assert_eq!(sequence.decode_int(), vec![1, 2, 3, 0]);
        assert_eq!(sequence.decode_dna(), b"CGTA");
    }

    #[test]
    fn inserts_base_at_beginning_middle_and_end() {
        let mut sequence = bitpacked::new("ACGT");

        sequence.insert(0, 4); // NACGT
        sequence.insert(2, 3); // NATCGT
        sequence.insert(sequence.len(), 0); // NATCGTA

        assert_eq!(sequence.len(), 7);
        assert_eq!(sequence.decode_int(), vec![4, 0, 3, 1, 2, 3, 0]);
        assert_eq!(sequence.decode_dna(), b"NATCGTA");
    }

    #[test]
    fn removes_base_at_beginning_middle_and_end() {
        let mut sequence = bitpacked::new("ACGTN");

        sequence.remove(0); // CGTN
        assert_eq!(sequence.decode_dna(), b"CGTN");

        sequence.remove(1); // CTN
        assert_eq!(sequence.decode_dna(), b"CTN");

        sequence.remove(sequence.len() - 1); // CT
        assert_eq!(sequence.decode_dna(), b"CT");

        assert_eq!(sequence.len(), 2);
    }

    #[test]
    fn iterates_over_encoded_bases() {
        let sequence = bitpacked::new("ACGTN");

        let bases: Vec<u8> = sequence.iter().collect();

        assert_eq!(bases, vec![0, 1, 2, 3, 4]);
    }

    #[test]
    fn iterates_in_reverse() {
        let sequence = bitpacked::new("ACGTN");

        let bases: Vec<u8> = sequence.iter().rev().collect();

        assert_eq!(bases, vec![4, 3, 2, 1, 0]);
    }
}