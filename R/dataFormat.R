### format: RDB, RDBES, CS (more options?)
## WGBIOPTIM_tool: ADV, SampleOptim, SDTool, BioSim (more options?)
### If RDBES: use as input RDBES_object and run before the function masterTable()

dataFormat <- function(format = NULL, input_data, WGBIOPTIM_tool)
{
  if (!is.null(format) && format=="RDBES")
    
  {
    ### ouput_data is a list of tables: sampleData, lengthData, biologicalData
    output_data <- masterTable(input_data)
    
  } else if (!is.null(format) && format=="RDB")
  
  {
    
    ### as 
    
  } else if (!is.null(format) && format=="CS")
    
  {
    ###
    
  } else { 
          output_data <- NULL
          message("Please define your input data format")
          }
    
  return(output_data)
}
