#with open("./COMPASS_PNNL/Compass_CRC_TR_302_WaterLevel600B.dat", 'r') as inputFile:
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
with open("/home/whf/E3SM/Teri/swamp-001.tec",'r') as inputFile:
    read_lines = inputFile.readlines()
    parameter_str = read_lines[1] #read in the 2nd line, getting variable names as strings
    parameter_list = parameter_str.replace('"','').replace('\n','').split(',')
    result_list = [] #define a list 
    dic = {} #define a dic conrresponding to result_list, each dic will have a list in result_list
    for i in range(0,len(parameter_list)): #loop for reading in data line by line
        dic = {parameter_list[i]:[]} #define a default format for each dic as variable_name:value_list
        result_list.append(dic) #append value list from each line in the data file
    for i in range(3, len(read_lines)): #start from the 4th line
        data_list = read_lines[i].split('  ') #begin reading data, split with two splaces '  '
        for j in range(0,len(result_list)): #loop for reading in data at each line
            for key,value in result_list[j].items(): #key is the variable name, value is the data
                value.append(data_list[j])            
    var_value1 = []
    var_value2 = []
    for i in range(0,len(result_list)):
       for key,value in result_list[i].items():
         if "Total NH4+" in key: #search for variable name and read in variable to array
            var_value1 = np.array(value, dtype=np.float32)
            print(var_value1)
         elif "Z [m]" in key: #read in depth to array
            var_value2 = np.array(value, dtype=np.float32)
            print(var_value2)
 #plotting variable in depth
    plt.figure(1);plt.clf()
    plt.plot(var_value1,var_value2)
    plt.xlabel('Total NH4+')
    plt.ylabel('Z [m]')
    ax = plt.gca()
    ax.invert_yaxis()
    #plt.ylim(max(var_value2), min(var_value2)) # this is another method to revert yaxis
    plt.show()    
