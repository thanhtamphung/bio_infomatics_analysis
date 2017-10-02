0/ This version is to convert and parse file vcf into table file so it is easy to handle, visualize, analyze 

1/ open file main.py, modify input file's path
    then
    run 
     $ python main.py

2/ or run directly without modification main.py

      python -c 'from main import main; main("path/of/your/vcf_file.vcf")'
  
  for example 
  
      python -c 'from main import main; main("sample_vcf_file.vcf")'
