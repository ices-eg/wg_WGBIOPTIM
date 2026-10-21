### format: RDB, RDBES, CS_CL (more options?)
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
    
  } else if (!is.null(format) && format=="CS_CL" && WGBIOPTIM_tool=="SDTool")
    
  {
       output_data1 <- SDTool::convert_cs_to_costCS(
      input_data[[1]],
      write_csv = TRUE,
      out_dir = getwd()
    )
    
    output_data2 <- SDTool::convert_cl_to_costCL(
      input_data[[2]],
      write_csv = TRUE,
      out_dir = getwd()
    )

    output_data=list(output_data1,output_data2)
    
  } else { 
          output_data <- NULL
          message("Please define your input data format")
          }
    
  return(output_data)
}
