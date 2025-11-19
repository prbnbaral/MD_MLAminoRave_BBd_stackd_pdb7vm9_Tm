def random_pick(x, x_dt, w1, w2, training_len):
    """ ramdomly pick (x, x_dt) pair from data set
    
        Parameters
        ----------
        x : np.array
            present trajectory.
         
        x_dt : np.array
            future trajectory.
            
        w1 : np.array
            reweighting factores in objective function before P(X_t | \chi )
            
        w2 : np.array
            reweighting factores in objective function before P(X | \chi )
            
        training_len: int
            length of the return data set
            
        
        Returns
        -------
        x1 : np.array
            ramdonly selected data pionts from present trajectory.
         
        x2 : np.array
            future trajectory corresponds to selected data points in x1.
            
        w1 : np.array
            coressponding reweighting factores in objective function before P(X_t | \chi )
        w1 : np.array
            coressponding reweighting factores in objective function before P(X | \chi )
    """
    indices = np.arange( np.shape(x)[0])
    np.random.shuffle(indices)
    indices = indices[:training_len]
    x = x[indices, :]
    x_dt = x_dt[indices, :]
    w1 = w1[indices]
    w2 = w2[indices]
    print('%i data points are used in this training'%len(indices))
    return x, x_dt, w1, w2
