def dynamic_correction_loss(x, w1, w2):
    """loss function with dynamic correction"""
    def custom_loss(y_true, y_pred ):
         ce1 = mean_squared_error(y_true, y_pred )
         ce2 = mean_squared_error(x, y_pred)  
         return (w1[:,0]*ce1+w2[:,0]*ce2)    
    return custom_loss
