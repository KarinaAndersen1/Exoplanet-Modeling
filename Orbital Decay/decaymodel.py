class DecayModel(object):
    def __init__(self, params, epoch, max_error = 1.0, nthreads = 1):
        self.epoch = epoch
        self.t0 = params.t0
        self.dp_dt = params.dp_dt
        self.P = params.P

    def calc_error(self):
        err = np.max(np.abs(t-t0))

    def decay_curve(self, params):
        self.t0 = params.t0
        self.dp_dt = params.dp_dt
        self.P = params.P

        t = self.t0 + self.P * self.epoch + 0.5 * self.P * self.dp_dt * self.epoch

        return t

class DecayParams(object):
    '''
    :param t0: Mid-transit time
    :type t0: Float

    :param E: Relative epoch
    :type E: Int

    :param P: Period
    :type P: Float

    :param dp_dt: Change in period (min)
    :type dp_dt: Float
    '''

    def __init__(self):
        self.t0 = None
        self.P = None
        self.dp_dt = None