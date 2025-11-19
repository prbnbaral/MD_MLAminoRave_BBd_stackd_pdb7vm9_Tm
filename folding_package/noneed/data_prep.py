def data_prep(system_name, number_trajs, predictive_step):
    """ Read the input trajectory files.
        Prepare x, x_t trajectory and corresponding reweighting factors
        
        Parameters
        ----------
        system_name : string
            Name of the sytem.
            
        number_trajs : int
            Number of trajectories.
        
        predictive_step : int
            Predictive time delay.
        
        Returns
        -------
        X : np.array
            present trajectory.
         
        Y : np.array
            future trajectory.
            
        W1 : np.array
            reweighting factores in objective function before P(X_t | \chi )
            
        W2 : np.array
            reweighting factores in objective function before P(X | \chi )
    """
    
    
    for j in range(number_trajs):
        traj_file_name = 'input/x_'+system_name+'_%i.npy'%j   #present trajecotry of the shape n*d, where n is the MD steps and d is the number of order parameters
        w_file_name = 'input/w_'+system_name+'_%i.npy'%j      #weights correspond to trajecotry in x. Calculated by exp(beta*V)   
        if predictive_step==0:
            x = np.load(traj_file_name)
            y = x[:,:]      
            w1 =  np.load(w_file_name)   
            w2 = np.zeros( np.shape(w1) )        
        else:
            x = np.load(traj_file_name)
            y = x[predictive_step: , :]
            x = x[:-predictive_step, :]
            w = np.load(w_file_name)
            w_x = w[:-predictive_step]
            w_y = w[predictive_step:]            
            w1 = ( w_x * w_y )**0.5
            w2 =  w_x**0.5*( w_x**0.5- w_y**0.5)
        try:
            X = np.append(X, x, axis = 0)
            Y = np.append(Y, y, axis = 0)
            W1 = np.append(W1, w1, axis = 0)
            W2 = np.append(W2, w2, axis = 0)
        except:
            X = x
            Y = y
            W1 = w1
            W2 = w2
    normaliztion_factor = np.sum(W1)/len(W1)  
    W1 /= normaliztion_factor
    W2 /= normaliztion_factor    
    print('length of data:%i'%np.shape(X)[0] )
    print('number of order parameters:%i'%np.shape(X)[1] )
    print('min reweighting factor:%f'%np.min(W1))
    print('max reweighting factor:%f'%np.max(W1))   
    return X, Y, W1, W2
