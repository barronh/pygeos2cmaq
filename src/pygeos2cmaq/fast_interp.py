import numpy as np

def get_interp_w(x, xn, axis = 0):
    """
    Returns array of weighting variables with shape (xn.shape, x.shape)
    """
    dx = np.ma.diff(x, axis = axis)
    if (dx > 0).all():
        reverse = False
    elif (dx < 0).all():
        reverse = True
    else:
        raise TypeError('Must be either ascending or descending')
    if reverse:
        idx = tuple([slice(None)] * axis + [slice(None, None, -1)])
        x = x[idx]
        xn = xn[idx]    
    dx = np.ma.diff(x, axis = axis)
    xtmp = np.ma.concatenate([x, np.max([xn.max(axis = axis, keepdims = True), x.max(axis = axis, keepdims = True) + dx.take([-1], axis = axis)], axis = 0)], axis = axis)
    dx = np.ma.diff(xtmp, axis = axis)
    
    xnshape = range(xn.ndim)
    xshape = range(x.ndim)
    xnshape.insert(axis + 1, -1)
    xshape.insert(axis, -1)
    dimiter = np.nditer([xn, x, dx, None], flags=['external_loop'], op_axes=[xnshape, xshape, xshape, None])
    for xn_, x_, dx_, pct_dx_ in dimiter:
        pct_dx_[...] = (xn_ - x_) / dx_
    pct_dx = dimiter.operands[3]
    #pct_dx = (xn[:, None] - x[None, :]) / dx[None, :]
    test = np.ma.logical_and(pct_dx < 1, pct_dx >= 0)
    #test[:, 1:] = np.ma.logical_or(test[:, 1:], test[:, :-1])
    mask = False == test
    w = np.ma.masked_where(mask, pct_dx)
    w = 1 - w
    idxs = tuple([slice(None)] * (axis + 1) + [slice(1, None)])
    idxe = tuple([slice(None)] * (axis + 1) + [slice(None, -1)])
    
    w[idxs] = np.ma.array([w[idxs], 1 - w[idxe]], dtype = w.dtype).max(0)
    nmax = xn.max(axis = axis, keepdims = True)
    omax = x.max(axis = axis, keepdims = True)
    if np.any(nmax > omax):
        # Needs special case extrapolation code
        idx = np.ma.where(xn >= x.max(axis = axis, keepdims = True))
        for thisidx in zip(*idx):
            idxs = list(thisidx)
            idxs.insert(axis + 1, -1)
            idxs = tuple(idxs)
            idxe = list(thisidx)
            idxe.insert(axis + 1, -2)
            idxe = tuple(idxe)
            w[idxs] = pct_dx[idxe]
            w[idxe] = 1 - w[idxs]
        
    
    nmin = xn.min(axis = axis, keepdims = True)
    omin = x.min(axis = axis, keepdims = True)
    if np.any(nmin < omin):
        # Needs special case extrapolation code
        idx = np.ma.where(xn <= x.min(axis = axis, keepdims = True))
        for thisidx in zip(*idx):
            idxs = list(thisidx)
            idxs.insert(axis + 1, 0)
            idxs = tuple(idxs)
            idxe = list(thisidx)
            idxe.insert(axis + 1, 1)
            idxe = tuple(idxe)
            w[idxs] = 1 - pct_dx[idxs]
            w[idxe] = 1 - w[idxs]
    try:
        np.testing.assert_allclose(w.sum(axis + 1).filled(1), 1, rtol=1e-6)
    except:
        print w.sum(axis + 1); print w
        raise
        
    if reverse:
        idx = tuple([slice(None)] * axis + [slice(None, None, -1), slice(None, None, -1)])
        w = w[idx]
    return w

if __name__ == '__main__':
    x = np.array([1.,2.,3.])
    xn = np.array([.5, 1.75, 5.])
    y = np.array([3.,6.,9.])
    yhat2 = (get_interp_w(x, xn) * y[None]).sum(1)
    yhat1 = np.array([  1.5 ,   5.25,  15.  ])
    assert((yhat1 == yhat2).all())

    xn = np.array([ 1. ,  4.5,  8. ])
    x = np.array([  0.,   1.,   2.,   3.,   4.,   5.,   6.,   7.,   8.,   9.,  10.])
    y = np.array([ 10.,  15.,  20.,  25.,  30.,  35.,  40.,  45.,  50.,  55.,  60.])
    w = get_interp_w(x, xn)
    yhat1 = np.interp(xn, x, y)
    yhat2 = (w * y).sum(1)
    assert((yhat1 == yhat2).all())
    sout = np.array([ 0.99800003,  0.99599999,  0.99400002,  0.99199998,  0.99000001,
        0.98799998,  0.986     ,  0.97799999,  0.96700001,  0.95499998,
        0.94099998,  0.92500001,  0.90700001,  0.88599998,  0.86299998,
        0.83700001,  0.80699998,  0.77399999,  0.736     ,  0.69300002,
        0.64399999,  0.58899999,  0.52600002,  0.456     ,  0.37599999,
        0.285     ,  0.183     ,  0.064     ,  0.        ], dtype='f')
    sin = np.array([  1.00000000e+00,   9.84843075e-01,   9.69642341e-01,
         9.54442620e-01,   9.39243913e-01,   9.24044251e-01,
         9.08845544e-01,   8.93647790e-01,   8.78450155e-01,
         8.63252401e-01,   8.48054707e-01,   8.32857966e-01,
         8.15130472e-01,   7.92337418e-01,   7.67011523e-01,
         7.41689622e-01,   7.16370702e-01,   6.91052794e-01,
         6.59407377e-01,   6.21441424e-01,   5.83480537e-01,
         5.45528471e-01,   5.07588446e-01,   4.69661266e-01,
         4.31751013e-01,   3.93863648e-01,   3.56002122e-01,
         3.12176079e-01,   2.65558690e-01,   2.25441828e-01,
         1.91474453e-01,   1.62709877e-01,   1.38300866e-01,
         1.17552295e-01,   9.99152735e-02,   8.49244073e-02,
         6.70357868e-02,   4.79747690e-02,   3.40429507e-02,
         2.39078291e-02,   1.44230574e-02,   6.60990505e-03,
         2.81022908e-03,   1.08988350e-03,   3.73901683e-04,
         1.00436344e-04,   0.00000000e+00], dtype='f')
    cin = np.array([  2.42377597e+19,   2.39810105e+19,   2.37205472e+19,
         2.34552328e+19,   2.31844605e+19,   2.29036233e+19,
         2.25581897e+19,   2.21745129e+19,   2.18441009e+19,
         2.15111402e+19,   2.11889547e+19,   2.08634882e+19,
         2.04968297e+19,   2.00062188e+19,   1.94680079e+19,
         1.89224082e+19,   1.83772165e+19,   1.78238751e+19,
         1.71083185e+19,   1.62871207e+19,   1.54073256e+19,
         1.45386718e+19,   1.36580488e+19,   1.28490006e+19,
         1.20072145e+19,   1.11628281e+19,   1.02943315e+19,
         9.27813761e+18,   8.16032791e+18,   7.22400855e+18,
         6.40136494e+18,   5.67638931e+18,   5.03833611e+18,
         4.37529789e+18,   3.74544733e+18,   3.20009038e+18,
         2.45259510e+18,   1.69040732e+18,   1.16418202e+18,
         7.99332515e+17,   4.48051194e+17,   1.94980864e+17,
         7.82326919e+16,   2.91829950e+16,   1.09473895e+16,
         3.96509448e+15,   9.48706602e+14], dtype='f')
    cout = np.array([  2.42038815e+19,   2.41700012e+19,   2.41361230e+19,
         2.41022449e+19,   2.40683667e+19,   2.40344864e+19,
         2.40006082e+19,   2.38637542e+19,   2.36744249e+19,
         2.34649613e+19,   2.32157460e+19,   2.29212814e+19,
         2.25115990e+19,   2.20082426e+19,   2.15057877e+19,
         2.09521968e+19,   2.03218248e+19,   1.96165233e+19,
         1.87998940e+19,   1.78664328e+19,   1.67750587e+19,
         1.55352461e+19,   1.40853982e+19,   1.25456564e+19,
         1.07530577e+19,   8.62649885e+18,   6.18777656e+18,
         2.33120407e+18,   9.48706602e+14], dtype='d')
    yhat1 = np.interp(sout[::-1], sin[::-1], cin[::-1])[::-1]
    w = get_interp_w(sin, sout)
    yhat2 = (w * cin[None, :]).sum(1)
    np.testing.assert_allclose(yhat1, cout, rtol=1e-07, atol=0, err_msg='', verbose=True)
    np.testing.assert_allclose(yhat1, yhat2, rtol=1e-07, atol=0, err_msg='', verbose=True)
    
    sin = np.array([ 0.99 ,  0.955,  0.885,  0.72 ,  0.45 ,  0.15 ])
    cin = np.array([  4.00000000e-03,   3.25700000e-03,   1.88300000e-03,
         7.02300000e-04,   3.17700000e-05,   3.97000000e-08])
    cout = np.array([  4.00000000e-03,   4.00000000e-03,   4.00000000e-03,
         4.00000000e-03,   4.00000000e-03,   3.95754234e-03,
         3.91508574e-03,   3.74525683e-03,   3.51174302e-03,
         3.25699967e-03,   2.98219970e-03,   2.66814309e-03,
         2.31482867e-03,   1.90262813e-03,   1.72557316e-03,
         1.53952373e-03,   1.32485078e-03,   1.08871083e-03,
         8.16792131e-04,   6.35247046e-04,   5.13558207e-04,
         3.76968745e-04,   2.20511835e-04,   4.66706673e-05,
         2.39431913e-05,   1.43183346e-05,   3.53003282e-06,
         3.97000000e-08,   3.97000000e-08])

    yhat1 = np.interp(sout[::-1], sin[::-1], cin[::-1])[::-1]
    w = get_interp_w(sin, sout)
    yhat2 = (w * cin[None, :]).sum(1)

    np.testing.assert_allclose(yhat1, cout, rtol=1e-07, atol=0, err_msg='All values should be exact as they were produced from this function and archived', verbose=True)
    np.testing.assert_allclose(yhat1[5:-2], yhat2[5:-2], rtol=1e-07, atol=0, err_msg='Within inputs, all values should be close', verbose=True)
    np.testing.assert_array_less(yhat1[:5], yhat2[:5], err_msg = 'Expecting extrapolation to cause higher values')
    np.testing.assert_array_less(yhat2[-2:], yhat1[-2:], err_msg = 'Expecting extrapolation to cause lower values')
    
    from scipy.interpolate import LinearNDInterpolator
    
    sin = np.array([[ 1 ,  0.955,  0.885,  0.72 ,  0.45 ,  0 ],
                    [ 1 ,  0.95,  0.89,  0.7 ,  0.48 ,  0 ]])
    cin = np.array([[  4.00000000e-03,   3.25700000e-03,   1.88300000e-03,
         7.02300000e-04,   3.17700000e-05,   3.97000000e-08],
                    [  4.00000000e-03,   3.25700000e-03,   1.88300000e-03,
         7.02300000e-04,   3.17700000e-05,   3.97000000e-08]])
    cout = np.array([[3.96697820e-03, 3.93395542e-03, 3.90093362e-03, 3.86791084e-03, 
      3.83488905e-03, 3.80186627e-03, 3.76884447e-03, 3.63675531e-03, 
      3.45513346e-03, 3.25699967e-03, 2.98219970e-03, 2.66814309e-03, 
      2.31482867e-03, 1.90262813e-03, 1.72557316e-03, 1.53952373e-03, 
      1.32485078e-03, 1.08871083e-03, 8.16792131e-04, 6.35247046e-04, 
      5.13558207e-04, 3.76968745e-04, 2.20511835e-04, 4.66706673e-05, 
      2.65521275e-05, 2.01355564e-05, 1.29433552e-05, 4.55245399e-06, 
      3.97000000e-08],
     [3.97028038e-03, 3.94055988e-03, 3.91084026e-03, 3.88111976e-03, 
      3.85140014e-03, 3.82167964e-03, 3.79196002e-03, 3.67307978e-03, 
      3.50962011e-03, 3.33129975e-03, 3.05089965e-03, 2.68450027e-03, 
      2.27230012e-03, 1.85814302e-03, 1.71521616e-03, 1.55364692e-03, 
      1.36722041e-03, 1.16215151e-03, 9.26011588e-04, 6.80965011e-04, 
      5.31619617e-04, 3.63987096e-04, 1.71971797e-04, 3.01834850e-05, 
      2.48951008e-05, 1.88795654e-05, 1.21368768e-05, 4.27040687e-06, 
      3.97000000e-08]])
    yhat1 =  np.array([np.interp(sout[::-1], sin[0, ::-1], cin[0, ::-1])[::-1],
                       np.interp(sout[::-1], sin[1, ::-1], cin[1, ::-1])[::-1]])
    w = get_interp_w(sin, sout[None, :].repeat(2, 0), axis = 1)
    yhat2 = (w * cin[:, None, :]).sum(2)
    np.testing.assert_allclose(yhat1, cout, rtol=1e-07, atol=0, err_msg='All values should be exact as they were produced from this function and archived', verbose=True)
    np.testing.assert_allclose(yhat1[:], yhat2[:], rtol=1e-07, atol=0, err_msg='Within inputs, all values should be close', verbose=True)
