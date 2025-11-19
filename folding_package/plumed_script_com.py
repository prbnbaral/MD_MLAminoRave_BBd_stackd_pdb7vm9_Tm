def plumed_script(OPs_r, RCs, iteration, T, bias_factor, label, sigma, height, pace, output_path):
    
    """Write the assigned parameter values into the params file near the corresponding parameter identicator"""
    
    try:
        file = open(output_path+'plumed_s_'+str(iteration)+'.dat',"w+")        
        file.seek(0)
        #length=5
        
        try:
            #file.write(section.rjust(2) + "  " + type.rjust(2) + "  " + parameter.rjust(2) + "  " + value.rjust(8) + "                    \n")
            #file.write("This script is automatically created. \n\n\n\n")
            file.write("""""")
            # Defines the OPs to be combined
            for i in range(0, len(OPs_r.columns)):
                           file.write(str(OPs_r.columns[i])+": DISTANCE ATOMS="
                                      +str(int(OPs_r.columns[i][OPs_r.columns[i].find('d')+len('d'):OPs_r.columns[i].rfind('-')])+1)+','
                                      +str(int(OPs_r.columns[i][OPs_r.columns[i].find('-')+len('-'):OPs_r.columns[i].rfind('')])+1)+'\n' 
                           )
            # Combines the defined OPs
            file.write("\nCOMBINE "+"LABEL=RC_1"+" ARG=")
            
            for i in range(0, len(OPs_r.columns)):
                file.write(str(OPs_r.columns[i]))
                if i!=len(OPs_r.columns)-1:
                    file.write(',')
                if i==len(OPs_r.columns)-1:
                    file.write(' COEFFICIENTS=')
            for i in range(1, RCs.shape[1]):
                file.write(str(RCs.iloc[9, i]))
                if i!=RCs.shape[1]-1:
                    file.write(',')
                if i==RCs.shape[1]-1:
                    file.write(' PERIODIC=NO\n')
                    file.write('\n\n\n\n')
            
            # Define the settings of metadynamics run

            file.write("METAD ")
            file.write("BIASFACTOR="+str(bias_factor)+" ")
            file.write("TEMP="+str(T)+" ")
            file.write("ARG="+label+" ")
            file.write("SIGMA="+str(sigma)+" ")
            file.write("HEIGHT="+str(height)+" ")
            file.write("PACE="+str(pace)+" ")
            file.write("LABEL=metadynamics ")
            file.write("\n\n\n\n")
            
            file.write("PRINT ARG=")
            for i in range(0, len(OPs_r.columns)):
                file.write(str(OPs_r.columns[i]))
                file.write(',')
            file.write(",metadynamics.bias ")
            file.write("STRIDE=1 ")
            file.write("FILE=BIASED_COLVAR\n\n\n\n")
            file.write("""""")


            
        finally:
            file.close()
    except IOError:
        pass
