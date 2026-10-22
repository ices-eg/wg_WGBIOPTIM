\name{dataFormat}
\alias{dataFormat}
\title{Convert input data for WGBIOPTIM tools}
\description{
Converts fisheries sampling data according to their input format and, for
CS/CL data, the selected WGBIOPTIM tool.
}
\usage{
dataFormat(format = NULL, input_data, WGBIOPTIM_tool)
}
\arguments{
  \item{format}{A single character string identifying the input format.
  Implemented options are \code{"RDBES"} and \code{"CS_CL"}.
  The \code{"RDB"} option is recognised but not yet implemented.
  The default is \code{NULL}. Values are case-sensitive.}
  \item{input_data}{Input data to be converted. For \code{"RDBES"}, an
  RDBES object accepted by \code{masterTable()}. For \code{"CS_CL"}, a list
  containing CS data as its first element and CL data as its second element.
  The second element is required for \code{"SDTool"} but is not used for
  \code{"BioSim"}. CS and CL inputs must meet the requirements of the
  corresponding conversion functions in the \pkg{SDTool} package.}
  \item{WGBIOPTIM_tool}{A single character string identifying the target tool.
  For \code{format = "CS_CL"}, implemented values are \code{"SDTool"} and
  \code{"BioSim"}. This argument has no default but is not evaluated in the
  \code{"RDBES"} branch. Other tool names, including \code{"ADV"} and
  \code{"SampleOptim"}, do not currently have dedicated conversion branches.}
}
\details{
Conversion depends on the selected format and tool:
\describe{
  \item{\code{format = "RDBES"}}{
  Calls \code{masterTable(input_data)} internally. Supply the original
  RDBES object; a separate call to \code{masterTable()} is not required
  before calling this function.
  }
  \item{\code{format = "CS_CL"}, \code{WGBIOPTIM_tool = "SDTool"}}{
  Converts \code{input_data[[1]]} with
  \code{SDTool::convert_cs_to_costCS()} and \code{input_data[[2]]} with
  \code{SDTool::convert_cl_to_costCL()}. Both calls use
  \code{write_csv = FALSE}, \code{write_rds = TRUE}, and
  \code{out_dir = getwd()}.
  }
  \item{\code{format = "CS_CL"}, \code{WGBIOPTIM_tool = "BioSim"}}{
  Converts only \code{input_data[[1]]} with
  \code{SDTool::convert_cs_to_costCS()}, using \code{write_csv = TRUE}
  and \code{out_dir = getwd()}. The \code{write_rds} argument is not
  supplied, so the converter's default applies.
  }
  \item{\code{format = "RDB"}}{
  The conversion branch is currently empty and does not assign an output.
  This option should not be used until its implementation is completed.
  }
}
For CS/CL conversions, the \pkg{SDTool} package must be installed.
For RDBES conversion, \code{masterTable()} must be available to this
function. Required columns and input validation are delegated to the
underlying conversion functions.

The function writes files in the current working directory when the
selected converter is called with file-writing options enabled.
It does not expose arguments to change the output directory or disable
these side effects.
}
\value{
The result depends on the selected branch:
\describe{
  \item{\code{"RDBES"}}{
  The object returned by \code{masterTable()}, expected to be a list
  containing \code{sampleData}, \code{lengthData}, and
  \code{biologicalData}.
  }
  \item{\code{"CS_CL"} with \code{"SDTool"}}{
  An unnamed list of length two. The first element is the result of
  \code{SDTool::convert_cs_to_costCS()} and the second is the result of
  \code{SDTool::convert_cl_to_costCL()}.
  }
  \item{\code{"CS_CL"} with \code{"BioSim"}}{
  The object returned by \code{SDTool::convert_cs_to_costCS()}.
  }
  \item{Unmatched input options}{
  \code{NULL}, with the message
  \dQuote{Please define your input data format}.
  This includes \code{format = NULL} and unsupported CS/CL target tools,
  provided the conditions can be evaluated.
  }
}
The \code{"RDB"} branch has no defined return value in the current
implementation and will normally raise an error because
\code{output_data} has not been created.
}
\note{
This documentation describes the current implementation.
Argument values are not explicitly validated. Use non-missing,
length-one character strings for \code{format} and
\code{WGBIOPTIM_tool}. Omitting \code{WGBIOPTIM_tool} when
\code{format = "CS_CL"} raises an error.
}
\examples{
\dontrun{
# Convert an RDBES object; masterTable() is called internally.
rdbes_tables <- dataFormat(
  format = "RDBES",
  input_data = RDBES_object,
  WGBIOPTIM_tool = "SDTool"
)

# Convert CS and CL data for SDTool; RDS files are written to getwd().
sdtool_data <- dataFormat(
  format = "CS_CL",
  input_data = list(cs_data, cl_data),
  WGBIOPTIM_tool = "SDTool"
)

# Convert CS data for BioSim; CSV output is written to getwd().
biosim_data <- dataFormat(
  format = "CS_CL",
  input_data = list(cs_data),
  WGBIOPTIM_tool = "BioSim"
)
}
}
\keyword{manip}

