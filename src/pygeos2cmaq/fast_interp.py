import numpy as np

def get_interp_w(x, xn):
    dx = np.diff(x)
    if (dx > 0).all():
        reverse = False
    elif (dx < 0).all():
        reverse = True
    else:
        raise TypeError('Must be either ascending or descending')
    if reverse:
        x = x[::-1]
        xn = xn[::-1]    
    
    dx = np.append(np.diff(x), xn.max() - x.max())
        
    pct_dx = (xn[:, None] - x[None, :]) / dx[None, :]
    test = np.logical_and(pct_dx < 1, pct_dx >= 0)
    #test[:, 1:] = np.logical_or(test[:, 1:], test[:, :-1])
    mask = False == test
    w = np.ma.masked_where(mask, pct_dx)
    w = 1 - w
    w[:, 1:] = np.ma.array([w[:, 1:], 1 - w[:, :-1]], dtype = w.dtype).max(0)
    nmax = xn.max()
    omax = x[-1]
    if nmax > omax:
        # Needs special case extrapolation code
        for i in np.where(xn >= x.max())[0]:
            w[i, -1] = pct_dx[i, -2]
            w[i, -2] = 1 - w[i, -1]
        
    
    nmin = xn.min()
    omin = x[0]
    if nmin < omin:
        # Needs special case extrapolation code
        for i in np.where(xn <= x.min())[0]:
            w[i, 0] = pct_dx[i, 0]
            w[i, 1] = 1 - w[i, 0]
    w = w.filled(0)
    if reverse:
        w = w[::-1, ::-1]
    return w

if __name__ == '__main__':
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

    np.testing.assert_allclose(yhat1, cout, rtol=1e-07, atol=0, err_msg='', verbose=True)
    np.testing.assert_allclose(yhat1, yhat2, rtol=1e-07, atol=0, err_msg='', verbose=True)
    