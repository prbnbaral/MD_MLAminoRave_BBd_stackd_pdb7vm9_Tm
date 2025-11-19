def scaling(x):
    """ make order parametes with mean 0 and variance 1
        return new order parameter and scaling factors
        
        Parameters
        ----------
        x : np.array
            order parameters
            
        Returns
        ----------
        x : np.array
            order parameters after rescaling
        
        std_x : np.array
            resclaing factors of each OPs
              
     """ 
   
    x = x-np.mean(x, axis =0)
    std_x = np.std(x, axis =0)
    return x/std_x, std_x
