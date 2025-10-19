## HPA XML Parser

HPA XML Parser is an integrated toolkit for parsing, analyzing, and exporting RNA expression data from the Human Protein Atlas (HPA) XML datasets. The project contains two complementary implementations:
	1.	A Shiny web application written in R, providing a graphical interface for uploading, parsing, exploring, and exporting HPA XML data.
	2.	A Python script for programmatic and large-scale XML parsing and conversion to Excel, suitable for data analysis pipelines and integration into Jupyter notebooks.

The primary goal is to streamline the extraction of rnaExpression data (for both human brain regions and tissue-level expression) and export the processed results into structured .xlsx files for downstream bioinformatics analysis.

⸻

### Overview

The HPA XML Parser takes as input an XML file from the [Human Protein Atlas](https://www.proteinatlas.org/) RNA expression datasets. These files typically contain rnaExpression elements describing expression levels (nTPM) across different tissues or brain regions.

The toolkit parses these complex XML structures into well-formatted tabular data with attributes such as sample ID, variable name, normalized TPM (nTPM), sex, age, region, and organ. It supports two assay types:
	•	humanBrain: RNA expression data across brain regions.
	•	tissue: RNA expression data across multiple tissue types.

⸻

### 1. Shiny Application (R)

#### Features

The R Shiny app provides an interactive dashboard with the following capabilities:
	•	XML Upload & Parsing: Upload an HPA XML file and parse its contents into structured tables.
	•	Data Exploration: View parsed humanBrain and tissue data directly in the browser with pagination and filtering.
	•	Download: Export processed data as an Excel file containing separate sheets for brain and tissue datasets.
	•	Theme Selector: Customize the UI theme dynamically.
	•	Progress Indicators: Loading spinners and modals for responsive interaction during parsing.

#### User Interface Structure

The dashboard is built with shinydashboard and includes five tabs:
	•	Upload & Parse: Upload XML file and start parsing.
	•	humanBrain Data: Explore brain-specific RNA expression data.
	•	Tissue Data: Explore tissue-level RNA expression data.
	•	Download: Export parsed data as an .xlsx file.
	•	Select Theme: Choose from a range of UI themes using shinythemes.

#### How to Run the App

    1.	Install required R packages:

```R
install.packages(c("shiny", "shinydashboard", "shinythemes", "XML", "openxlsx", "DT", "shinycssloaders"))
```

    2.	Run the Shiny application:

```R
shiny::runApp("app_directory")
```

Replace "app_directory" with the path where the app.R file is located.

    3.	The application will open in your browser.

	4.	Follow these steps in the UI:
	    •	Upload the .xml file.
	    •	Click “Parse XML” to process the file.
	    •	Explore parsed data in the humanBrain and Tissue tabs.
	    •	Download the results as an Excel workbook.

⸻

### 2. Python Script

#### Features

The Python script provides a non-interactive way to parse the same HPA XML data for automated pipelines or integration with analysis workflows. It:

	•	Iterates through rnaExpression nodes.
	•	Extracts data for all regions or all tissues.
	•	Creates a comprehensive DataFrame for each assay type.
	•	Exports the combined dataset to an .xlsx file.

#### How It Works

The script uses xml.etree.ElementTree to traverse the XML structure and pandas for data handling and export. It defines two main functions:

	•	hpa_parser_all_all_fast(regionID, rt) – Parses all RNA samples for a specific region or tissue.
	•	loop_parser(path, file) – Iterates over all regions or tissues, aggregates data, and writes to Excel.

#### Example Usage

	1.	Install the required Python packages:

```bash
pip install pandas openpyxl
```

    2.	Run the script:

```python
python HPA_XML_Parser.py
```

	3.	Enter the input XML file path when prompted:

```plaintext
Insira o caminho do arquivo aqui: /path/to/your/file.xml
```

	4.	The script will automatically create a corresponding .xlsx file in the same directory.

#### Output Format

The resulting Excel file contains columns such as:

	•	ID: Sample identifier
	•	Variable: Variable name
	•	nTPM: Normalized expression level
	•	Sex: Donor sex
	•	Age: Donor age
	•	Region: Brain region or tissue name
	•	Organ: Organ name

⸻

#### File Structure

```
.
├── app.R                      # Shiny web application
├── HPA_XML_Parser.py         # Python script for batch parsing
├── README.md                 # Project documentation
└── example_data/            # Optional directory with example XML files
```

⸻

## Use Cases

This project is designed for researchers working with transcriptomic data from the Human Protein Atlas. Typical use cases include:

	•	Rapid exploration of RNA expression patterns in brain regions or tissues.
	•	Integration into data pipelines for downstream bioinformatics or machine learning analyses.
	•	Automated batch processing of multiple XML datasets.
	•	Generating structured data for visualization, statistical analysis, or cross-species comparisons.

⸻

## Requirements

### R Environment
	•	R ≥ 4.0
	•	shiny, shinydashboard, shinythemes, XML, openxlsx, DT, shinycssloaders

### Python Environment
	•	Python ≥ 3.8
	•	pandas, openpyxl

⸻

## License

This project is distributed under the MIT License. See the LICENSE file for details.

⸻

## Citation

If you use this tool in your research, please cite it as:

HPA XML Parser: A toolkit for structured extraction of RNA expression data from Human Protein Atlas XML datasets.
[GitHub Repository](https://github.com/drhrf/XML_HPA_parser) 