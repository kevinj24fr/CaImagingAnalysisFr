pkgname <- "CaImagingAnalysisFr"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('CaImagingAnalysisFr')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("calcium_correction")
### * calcium_correction

flush(stderr()); flush(stdout())

### Name: calcium_correction
### Title: Calcium Trace Correction
### Aliases: calcium_correction

### ** Examples

# Basic usage
raw <- generate_synthetic_data(3, 500)
corrected <- calcium_correction(raw)

# Custom parameters
corrected <- calcium_correction(raw, span = 0.3, normalize = FALSE)

# Legacy method
corrected <- calcium_correction(raw, method = "legacy")




cleanEx()
nameEx("calcium_pipeline")
### * calcium_pipeline

flush(stderr()); flush(stdout())

### Name: calcium_pipeline
### Title: Calcium Imaging Analysis Pipeline using targets
### Aliases: calcium_pipeline

### ** Examples

# Basic pipeline
pipeline <- calcium_pipeline()

# Custom parameters
pipeline <- calcium_pipeline(
  raw_data_path = "path_to_custom_data.csv",
  config = list(
    correction_method = "legacy",
    spike_method = "caiman"
  )
)




cleanEx()
nameEx("generate_synthetic_data")
### * generate_synthetic_data

flush(stderr()); flush(stdout())

### Name: generate_synthetic_data
### Title: Generate Synthetic Calcium Imaging Dataset
### Aliases: generate_synthetic_data

### ** Examples

# Basic usage
df <- generate_synthetic_data(3, 500)
head(df)

# Custom parameters
df <- generate_synthetic_data(n_cells = 10, n_time = 2000, spike_prob = 0.01)




cleanEx()
nameEx("infer_spikes")
### * infer_spikes

flush(stderr()); flush(stdout())

### Name: infer_spikes
### Title: Spike Inference Functions
### Aliases: infer_spikes

### ** Examples

# Basic usage
raw <- generate_synthetic_data(1, 500)
corrected <- calcium_correction(raw)
spikes <- infer_spikes(corrected$Cell_1)

# Try deep learning method (requires model)
# spikes <- infer_spikes(corrected$Cell_1, method = "deep", model_path = "cascade_model.h5")




cleanEx()
nameEx("launch_interactive_viewer")
### * launch_interactive_viewer

flush(stderr()); flush(stdout())

### Name: launch_interactive_viewer
### Title: Launch Interactive Viewer for Calcium Traces
### Aliases: launch_interactive_viewer

### ** Examples

# Basic usage
raw <- generate_synthetic_data(3, 500)
launch_interactive_viewer(raw)

# Custom port
launch_interactive_viewer(raw, port = 8080)




cleanEx()
nameEx("plot_cell_trace")
### * plot_cell_trace

flush(stderr()); flush(stdout())

### Name: plot_cell_trace
### Title: Plot Corrected Trace with Detected Spikes
### Aliases: plot_cell_trace

### ** Examples

# Basic usage
raw <- generate_synthetic_data(3, 500)
corrected <- calcium_correction(raw)
p <- plot_cell_trace(corrected, "Cell_1")
print(p)

# Custom options
p <- plot_cell_trace(corrected, "Cell_1", method = "caiman", show_spikes = FALSE)




cleanEx()
nameEx("plot_multiple_cells")
### * plot_multiple_cells

flush(stderr()); flush(stdout())

### Name: plot_multiple_cells
### Title: Plot Multiple Cell Traces
### Aliases: plot_multiple_cells

### ** Examples

raw <- generate_synthetic_data(6, 500)
corrected <- calcium_correction(raw)
p <- plot_multiple_cells(corrected, ncol = 2)




### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
