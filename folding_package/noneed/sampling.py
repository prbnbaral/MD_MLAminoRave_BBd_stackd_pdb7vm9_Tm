def sampling(args):
    """Sample the latent variable from a Normal distribution."""
    s_mean= args
    epsilon = K.random_normal(shape=(batch_size,rc_dim), mean=0.0, stddev=s_vari )
    s_noise = s_mean +  epsilon
    return s_noise
