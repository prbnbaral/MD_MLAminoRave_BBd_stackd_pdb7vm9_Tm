class WeightsHistory(Callback):
    def on_train_begin(self, logs={}):
        self.losses = []
        self.losses_vali = []
        self.weights0 = []
    def on_epoch_end(self, epoch, logs={}):
        self.losses.append(logs.get('loss'))
        self.losses_vali.append(logs.get('val_loss'))
        self.weights0.append( prave.layers[1].get_weights())
