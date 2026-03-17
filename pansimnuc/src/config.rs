use std::collections::HashMap;
use std::fs::File;
use std::io::{self, BufRead, BufReader};

#[derive(Debug, Clone)]
pub struct Config {
    pub sections: HashMap<String, HashMap<String, String>>,
}

impl Config {
    /// Parse a config file with sections marked by [header]
    pub fn from_file(path: &str) -> io::Result<Self> {
        let file = File::open(path)?;
        let reader = BufReader::new(file);
        let mut sections: HashMap<String, HashMap<String, String>> = HashMap::new();
        let mut current_section = "default".to_string();

        for line in reader.lines() {
            let line = line?;
            let trimmed = line.trim();

            // Skip empty lines and comments
            if trimmed.is_empty() || trimmed.starts_with('#') {
                continue;
            }

            // Check for section header
            if trimmed.starts_with('[') && trimmed.ends_with(']') {
                current_section = trimmed[1..trimmed.len() - 1].to_string();
                sections.entry(current_section.clone()).or_insert_with(HashMap::new);
                continue;
            }

            // Parse key=value pairs
            if let Some(eq_pos) = trimmed.find('=') {
                let key = trimmed[..eq_pos].trim().to_string();
                let value = trimmed[eq_pos + 1..].trim().to_string();
                sections
                    .entry(current_section.clone())
                    .or_insert_with(HashMap::new)
                    .insert(key, value);
            }
        }

        Ok(Config { sections })
    }

    /// Get a value from a specific section
    pub fn get(&self, section: &str, key: &str) -> Option<String> {
        self.sections
            .get(section)
            .and_then(|section_map| section_map.get(key).cloned())
    }

    /// Get all keys in a section
    pub fn keys_in_section(&self, section: &str) -> Vec<String> {
        self.sections
            .get(section)
            .map(|section_map| section_map.keys().cloned().collect())
            .unwrap_or_default()
    }

    /// Check if a section exists
    pub fn has_section(&self, section: &str) -> bool {
        self.sections.contains_key(section)
    }

    /// Flatten all config values into a single HashMap with keys like "section.key"
    pub fn flatten(&self) -> HashMap<String, String> {
        let mut flat = HashMap::new();
        for (section, values) in &self.sections {
            for (key, value) in values {
                flat.insert(format!("{}.{}", section, key), value.clone());
            }
        }
        flat
    }

    /// Get all values from all sections as a flat HashMap
    pub fn to_flat_hashmap(&self) -> HashMap<String, String> {
        self.flatten()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    fn create_test_config(content: &str) -> String {
        let temp_dir = std::env::temp_dir();
        let temp_path = temp_dir.join(format!(
            "test_config_{}_{}.conf",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        let mut file = File::create(&temp_path).expect("Failed to create test config file");
        file.write_all(content.as_bytes())
            .expect("Failed to write test config file");
        drop(file);
        temp_path.to_string_lossy().to_string()
    }

    #[test]
    fn test_config_parse_simple() {
        let content = "[database]
host=localhost
port=5432

[output]
format=json";

        let path = create_test_config(content);
        let config = Config::from_file(&path);
        assert!(config.is_ok());

        let config = config.unwrap();
        assert!(config.has_section("database"));
        assert!(config.has_section("output"));

        assert_eq!(config.get("database", "host"), Some("localhost".to_string()));
        assert_eq!(config.get("database", "port"), Some("5432".to_string()));
        assert_eq!(config.get("output", "format"), Some("json".to_string()));

        let _ = std::fs::remove_file(&path);
    }

    #[test]
    fn test_config_with_comments() {
        let content = "# This is a comment
[settings]
# Another comment
option1=value1
option2=value2";

        let path = create_test_config(content);
        let config = Config::from_file(&path);
        assert!(config.is_ok());

        let config = config.unwrap();
        assert_eq!(config.get("settings", "option1"), Some("value1".to_string()));
        assert_eq!(config.get("settings", "option2"), Some("value2".to_string()));

        let _ = std::fs::remove_file(&path);
    }

    #[test]
    fn test_config_input_with_earlgrey_gff_file() {
        let content = "[input]
gff_file=/tmp/main.gff3
fasta_file=/tmp/ref.fa
earlgrey_gff_file=/tmp/te.gff";

        let path = create_test_config(content);
        let config = Config::from_file(&path);
        assert!(config.is_ok());

        let config = config.unwrap();
        assert_eq!(config.get("input", "gff_file"), Some("/tmp/main.gff3".to_string()));
        assert_eq!(config.get("input", "fasta_file"), Some("/tmp/ref.fa".to_string()));
        assert_eq!(config.get("input", "earlgrey_gff_file"), Some("/tmp/te.gff".to_string()));

        let _ = std::fs::remove_file(&path);
    }

    #[test]
    fn test_config_keys_in_section() {
        let content = "[features]
exon=true
intron=false
intergenic=true";

        let path = create_test_config(content);
        let config = Config::from_file(&path);
        assert!(config.is_ok());

        let config = config.unwrap();
        let keys = config.keys_in_section("features");
        assert_eq!(keys.len(), 3);
        assert!(keys.contains(&"exon".to_string()));
        assert!(keys.contains(&"intron".to_string()));
        assert!(keys.contains(&"intergenic".to_string()));

        let _ = std::fs::remove_file(&path);
    }

    #[test]
    fn test_config_flatten_with_rates() {
        let content = "[exons]
mutation_rate=1e-8
selection_coefficient=0.1

[introns]
mutation_rate=1e-9
selection_coefficient=0.05

[intergenic]
mutation_rate=5e-9
selection_coefficient=0.02";

        let path = create_test_config(content);
        let config = Config::from_file(&path);
        assert!(config.is_ok());

        let config = config.unwrap();
        let flat = config.flatten();

        // Check that all flattened keys exist
        assert_eq!(flat.get("exons.mutation_rate"), Some(&"1e-8".to_string()));
        assert_eq!(flat.get("exons.selection_coefficient"), Some(&"0.1".to_string()));
        assert_eq!(flat.get("introns.mutation_rate"), Some(&"1e-9".to_string()));
        assert_eq!(flat.get("introns.selection_coefficient"), Some(&"0.05".to_string()));
        assert_eq!(flat.get("intergenic.mutation_rate"), Some(&"5e-9".to_string()));
        assert_eq!(flat.get("intergenic.selection_coefficient"), Some(&"0.02".to_string()));

        // Total of 6 values
        assert_eq!(flat.len(), 6);

        let _ = std::fs::remove_file(&path);
    }
}
