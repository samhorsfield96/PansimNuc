use std::collections::HashMap;
use std::fs::File;
use std::io::{self, BufRead, BufReader};

#[derive(Debug, Clone)]
pub struct Config {
    pub sections: HashMap<String, HashMap<String, String>>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct PopulationSplitConfig {
    pub population_splits: Vec<usize>,
    pub generation_splits: Vec<usize>,
    pub migration_rate: f64,
}

impl PopulationSplitConfig {
    pub fn new() -> Self {
        PopulationSplitConfig {
            population_splits: Vec::new(),
            generation_splits: Vec::new(),
            migration_rate: 0.0,
        }
    }
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
                sections
                    .entry(current_section.clone())
                    .or_insert_with(HashMap::new);
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

    /// Parse a comma-separated list of positive integers from a section key.
    pub fn get_usize_vec(&self, section: &str, key: &str) -> Result<Vec<usize>, String> {
        let value = self
            .get(section, key)
            .ok_or_else(|| format!("Missing required config key: {}.{}", section, key))?;

        if value.trim().is_empty() {
            return Err(format!(
                "Config key '{}.{}' must contain at least one integer",
                section, key
            ));
        }

        value
            .split(',')
            .map(|part| {
                let trimmed = part.trim();
                if trimmed.is_empty() {
                    return Err(format!(
                        "Config key '{}.{}' contains an empty list entry",
                        section, key
                    ));
                }
                trimmed.parse::<usize>().map_err(|_| {
                    format!(
                        "Config key '{}.{}' contains a non-integer value: '{}'",
                        section, key, trimmed
                    )
                })
            })
            .collect()
    }

    /// Parse an f64 from a section key.
    pub fn get_f64(&self, section: &str, key: &str) -> Result<f64, String> {
        let value = self
            .get(section, key)
            .ok_or_else(|| format!("Missing required config key: {}.{}", section, key))?;

        value.parse::<f64>().map_err(|_| {
            format!(
                "Config key '{}.{}' must be a floating-point number",
                section, key
            )
        })
    }

    /// Parse and validate migration split settings from [population].
    pub fn population_split_config(&self) -> Result<PopulationSplitConfig, String> {
        let population_splits = self.get_usize_vec("population", "population_splits")?;
        let generation_splits = self.get_usize_vec("population", "generation_splits")?;

        if population_splits.len() != generation_splits.len() {
            return Err(format!(
                "population.population_splits and population.generation_splits must have the same length (got {} and {})",
                population_splits.len(),
                generation_splits.len()
            ));
        }

        let migration_rate = self.get_f64("population", "migration_rate")?;

        Ok(PopulationSplitConfig {
            population_splits,
            generation_splits,
            migration_rate,
        })
    }

    /// Parse tracking regions from [tracking] into a list of (contig, start, end) tuples.
    pub fn tracking_regions(&self) -> Result<Vec<(usize, usize, usize)>, String> {
        let contigs: Vec<usize> = self
            .get_usize_vec("tracking", "contig")?
            .into_iter()
            .collect();

        let starts: Vec<usize> = self.get_usize_vec("tracking", "start")?;
        let ends: Vec<usize> = self.get_usize_vec("tracking", "end")?;

        if contigs.len() != starts.len() || contigs.len() != ends.len() {
            return Err(format!(
                "tracking.contig, tracking.start, and tracking.end must all have the same number of entries (got {}, {}, {})",
                contigs.len(), starts.len(), ends.len()
            ));
        }

        Ok(contigs.into_iter().zip(starts).zip(ends).map(|((c, s), e)| (c, s, e)).collect())
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
    fn test_tracking_regions_parses_correctly() {
        let content = "[tracking]
contig=0,0
start=0,60000
end=20000,100000
";
        let path = create_test_config(content);
        let config = Config::from_file(&path).unwrap();
        let regions = config.tracking_regions().unwrap();
        assert_eq!(regions, vec![
            (0, 0, 20000),
            (0, 60000, 100000),
        ]);
    }

    #[test]
    fn test_tracking_regions_mismatched_lengths_returns_error() {
        let content = "[tracking]
contig=0,0,1
start=0,60000
end=20000,100000
";
        let path = create_test_config(content);
        let config = Config::from_file(&path).unwrap();
        assert!(config.tracking_regions().is_err());
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

        assert_eq!(
            config.get("database", "host"),
            Some("localhost".to_string())
        );
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
        assert_eq!(
            config.get("settings", "option1"),
            Some("value1".to_string())
        );
        assert_eq!(
            config.get("settings", "option2"),
            Some("value2".to_string())
        );

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
        assert_eq!(
            config.get("input", "gff_file"),
            Some("/tmp/main.gff3".to_string())
        );
        assert_eq!(
            config.get("input", "fasta_file"),
            Some("/tmp/ref.fa".to_string())
        );
        assert_eq!(
            config.get("input", "earlgrey_gff_file"),
            Some("/tmp/te.gff".to_string())
        );

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
        assert_eq!(
            flat.get("exons.selection_coefficient"),
            Some(&"0.1".to_string())
        );
        assert_eq!(flat.get("introns.mutation_rate"), Some(&"1e-9".to_string()));
        assert_eq!(
            flat.get("introns.selection_coefficient"),
            Some(&"0.05".to_string())
        );
        assert_eq!(
            flat.get("intergenic.mutation_rate"),
            Some(&"5e-9".to_string())
        );
        assert_eq!(
            flat.get("intergenic.selection_coefficient"),
            Some(&"0.02".to_string())
        );

        // Total of 6 values
        assert_eq!(flat.len(), 6);

        let _ = std::fs::remove_file(&path);
    }

    #[test]
    fn test_population_split_config_parses_successfully() {
        let content = "[population]
population_splits=1,2,1,3
generation_splits=1,5,10,15
migration_rate=0.01";

        let path = create_test_config(content);
        let config = Config::from_file(&path).unwrap();
        let split_config = config.population_split_config().unwrap();

        assert_eq!(split_config.population_splits, vec![1, 2, 1, 3]);
        assert_eq!(split_config.generation_splits, vec![1, 5, 10, 15]);
        assert_eq!(split_config.migration_rate, 0.01);

        let _ = std::fs::remove_file(&path);
    }

    #[test]
    fn test_population_split_config_requires_equal_lengths() {
        let content = "[population]
population_splits=1,2,1
generation_splits=1,5
migration_rate=0.01";

        let path = create_test_config(content);
        let config = Config::from_file(&path).unwrap();
        let result = config.population_split_config();

        assert!(result.is_err());
        let err = result.err().unwrap();
        assert!(err.contains("must have the same length"));

        let _ = std::fs::remove_file(&path);
    }

    #[test]
    fn test_population_split_config_rejects_non_integer_split_values() {
        let content = "[population]
population_splits=1,2,a
generation_splits=1,5,10
migration_rate=0.01";

        let path = create_test_config(content);
        let config = Config::from_file(&path).unwrap();
        let result = config.population_split_config();

        assert!(result.is_err());
        let err = result.err().unwrap();
        assert!(err.contains("non-integer"));

        let _ = std::fs::remove_file(&path);
    }
}
