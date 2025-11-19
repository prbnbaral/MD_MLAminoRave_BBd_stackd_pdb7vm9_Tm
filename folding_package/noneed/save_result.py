def save_result(system_name, op_dim, time_delay, trials, s_vari, training_size, network_info, save_path):
    ''' save final result (linear combinaton coefficients of OPs) to a txt file
        Parameters
        ----------
        system_name : string
            Name of the system.
            
        op_dim : int
            Dimensionality of order parameters.
        
        time_delay : int
            Predictive time delay.
        
        trials: list
            Indexes of all trials.
            
        s_vari: float
            Variance of noise added to the decoder.
        
        training_size: int
            Total number of training data points.
        
        network_info: string
            Other detail of neural network.
        
        save_path: string
            Directory of where the final result is saved.
            
        Returns
        ----------
            None  
            Result is saved to a txt file.
    
    '''
    weights = []
    for dt in time_delay:
        Loss = []
        Weights = []
        for trial in trials: 
            save_dir = system_name+'_dt'+str(dt)+'_trail'+str(trial)+'_svar'+str(s_vari)+'_train_size'+str(training_size)+network_info+'.npy'
            Result_loss = np.load(save_path+'Loss_'+save_dir) 
            Result_weights = np.load(save_path+'Weights_'+save_dir) 
            Loss.append(np.average( Result_loss[-2:,-1] ))
            Weights.append( Result_weights[-1,:,:] )
        
        Weights = np.array( Weights )
        min_index = np.argmin(Loss)
        weights.append( Weights[min_index,:,:]  )
    weights = np.array(weights)
    
    ###save weights vs. time delay###
    head = 'time_delay/MD_step  '
    nunmber_rcs = np.shape(Weights)[-1]
    print('There are %i reaction coordinates'%nunmber_rcs)
    for j in range(op_dim):
        head+='op%i '%(j+1)
    for j in range(len(time_delay)):
         result_given_dt = np.concatenate((np.transpose( [[time_delay[j]]*nunmber_rcs] ), np.transpose(weights[j,:,:])), axis =-1)
         try:
             final_result = np.concatenate((final_result, result_given_dt), axis=0)
         except:
             final_result = result_given_dt
            
    np.savetxt( save_path+'final_result_'+system_name+'_svar'+str(s_vari)+'_train_size'+str(training_size)+network_info+'.txt', final_result, header=head, comments='###', newline='\n')


