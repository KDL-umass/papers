from dragonnet.src.experiment.models import *
import os
import glob
import argparse
import tensorflow as tf
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import train_test_split
import keras.backend as K
from keras.optimizers import rmsprop, SGD, Adam
from keras.callbacks import EarlyStopping, ModelCheckpoint, TensorBoard, ReduceLROnPlateau, TerminateOnNaN
import warnings
from dragonnet.src.semi_parametric_estimation.ate import *
import numpy as np
tf.compat.v1.logging.set_verbosity(tf.compat.v1.logging.ERROR)
warnings.filterwarnings("ignore", message=r"Passing", category=FutureWarning)
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '2'


def sigmoid(x):
  return 1 / (1 + np.exp(-x))

def _split_output(yt_hat, t, y,  x, binary_flag = False):
    q_t0 = yt_hat[:, 0].copy()
    q_t1 = yt_hat[:, 1].copy()
    
    g = yt_hat[:, 2].copy()
    if binary_flag:
        q_t0 = sigmoid(q_t0)
        q_t1 = sigmoid(q_t1)
    
    if yt_hat.shape[1] == 4:
        eps = yt_hat[:, 3][0]
    else:
        eps = np.zeros_like(yt_hat[:, 2])

    y = y.copy()
    var = "average propensity for treated: {} and untreated: {}".format(g[t.squeeze() == 1.].mean(),
                                                                        g[t.squeeze() == 0.].mean())
    print(var)
    
    y_pred = np.zeros_like(y)
    for i in range(t.shape[0]):
        if t[i] == 1.:
            y_pred[i] = q_t1[i]
        else:
            y_pred[i] = q_t0[i]

    return {'q_t0': q_t0, 'q_t1': q_t1, 'g': g, 't': t, 'y': y, 'x': x, 'eps': eps, 'y_pred': y_pred}


def train_and_predict_dragons(train_outcome, train_treatment, train_covariates, test_outcome, test_treatment, test_covariates, 
                              targeted_regularization=True, knob_loss=dragonnet_loss_binarycross, ratio=1., dragon='', val_split=0.2, batch_size=64, binary_flag = False):
    verbose = 0
    # y_scaler = StandardScaler().fit(train_outcome)
    # y_train = y_scaler.transform(train_outcome)
    # y_test = y_scaler.transform(test_outcome)

    train_outputs = []
    test_outputs = []
   
    if dragon == 'tarnet':
        dragonnet = make_tarnet(train_covariates.shape[1], 0.01, binary_flag = binary_flag)

    elif dragon == 'dragonnet':
        print("I am here making dragonnet")
        dragonnet = make_dragonnet(train_covariates.shape[1], 0.01, binary_flag = binary_flag)

    metrics = [regression_loss, binary_classification_loss, treatment_accuracy, track_epsilon]

    if targeted_regularization:
        loss = make_tarreg_loss(ratio=ratio, dragonnet_loss=knob_loss, binary_flag = binary_flag)
    else:
        loss = knob_loss(binary_flag = binary_flag)

    # for reporducing the IHDP experimemt

    i = 0
    tf.random.set_random_seed(i)
    np.random.seed(i)
    x_train, x_test = train_covariates, test_covariates
    t_train, t_test = train_treatment, test_treatment
    y_train, y_test = train_outcome, test_outcome
    yt_train = np.concatenate([y_train, t_train], 1)
    
    import time;
    start_time = time.time()

    dragonnet.compile(
        optimizer=Adam(lr=1e-3),
        loss=loss, metrics=metrics)

    adam_callbacks = [
        TerminateOnNaN(),
        EarlyStopping(monitor='val_loss', patience=2, min_delta=0.),
        ReduceLROnPlateau(monitor='loss', factor=0.5, patience=5, verbose=verbose, mode='auto',
                          min_delta=1e-8, cooldown=0, min_lr=0)

    ]

    dragonnet.fit(x_train, yt_train, callbacks=adam_callbacks,
                  validation_split=val_split,
                  epochs=100,
                  batch_size=batch_size, verbose=verbose)

    # sgd_callbacks = [
    #     TerminateOnNaN(),
    #     EarlyStopping(monitor='val_loss', patience=40, min_delta=0.),
    #     ReduceLROnPlateau(monitor='loss', factor=0.5, patience=5, verbose=verbose, mode='auto',
    #                       min_delta=0., cooldown=0, min_lr=0)
    # ]

    # sgd_lr = 1e-5
    # momentum = 0.9
    # dragonnet.compile(optimizer=SGD(lr=sgd_lr, momentum=momentum, nesterov=True), loss=loss,
    #                   metrics=metrics)
    # dragonnet.fit(x_train, yt_train, callbacks=sgd_callbacks,
    #               validation_split=val_split,
    #               epochs=300,
    #               batch_size=batch_size, verbose=verbose)

    # elapsed_time = time.time() - start_time
    # print("***************************** elapsed_time is: ", elapsed_time)

    yt_hat_test = dragonnet.predict(x_test)
    yt_hat_train = dragonnet.predict(x_train)

    test_outputs += [_split_output(yt_hat_test, t_test, y_test, x_test, binary_flag = binary_flag)]
    train_outputs += [_split_output(yt_hat_train, t_train, y_train, x_train, binary_flag = binary_flag)]
    K.clear_session()

    return test_outputs, train_outputs


def run_data(train_outcome, train_treatment, train_covariates, test_outcome, test_treatment, test_covariates, 
             knob_loss=dragonnet_loss_binarycross,
             ratio=1., dragon='', binary_flag = False):
    print("the dragon is {}".format(dragon))
    final_result = {}
    for is_targeted_regularization in [True]:
        if is_targeted_regularization:
            key = "t_reg_true"
        else:
            key = "t_reg_false"
        final_result[key] = {}
        results = {}
        print("Is targeted regularization: {}".format(is_targeted_regularization))
        test_outputs, train_output = train_and_predict_dragons(train_outcome, train_treatment, train_covariates, test_outcome, test_treatment, test_covariates, 
                                                                    targeted_regularization=is_targeted_regularization,
                                                                    knob_loss=knob_loss, ratio=ratio, dragon=dragon,
                                                                    val_split=0.2, batch_size=64, binary_flag = binary_flag)

        train_results, test_results = train_output[0], test_outputs[0]
        
        qt0_tr, qt1_tr, g_tr, t_tr, y_tr = train_results['q_t0'].reshape(-1, 1), train_results['q_t1'].reshape(-1, 1), train_results['g'].reshape(-1, 1), train_results['t'].reshape(-1, 1), train_results['y_pred'].reshape(-1, 1)
        psi_n_tr = psi_naive(qt0_tr, qt1_tr, g_tr, t_tr, y_tr, truncate_level=0.0)
        # psi_tmle_tr, psi_tmle_std_tr, eps_hat_tr, initial_loss_tr, final_loss_tr, g_loss_tr = psi_tmle_cont_outcome(qt0_tr, qt1_tr, g_tr, t_tr,
        #                                                                                       y_tr,
        #                                                                                       truncate_level=0.0)    
            
        qt0_te, qt1_te, g_te, t_te, y_te = test_results['q_t0'].reshape(-1, 1), test_results['q_t1'].reshape(-1, 1), test_results['g'].reshape(-1, 1), test_results['t'].reshape(-1, 1), test_results['y_pred'].reshape(-1, 1)
        psi_n_te = psi_naive(qt0_te, qt1_te, g_te, t_te, y_te, truncate_level=0.0)
        # psi_tmle_te, psi_tmle_std_te, eps_hat_te, initial_loss_te, final_loss_te, g_loss_te = psi_tmle_cont_outcome(qt0_te, qt1_te, g_te, t_te,
        #                                                                                       y_te,
        #                                                                                       truncate_level=0.0)                                    
        results["train_ate"] = psi_n_tr
        results["train_outcome"] = y_tr
        results["test_ate"] = psi_n_te
        results["test_outcome"] = y_te
       
        final_result[key] = results
    return final_result


def run(train_outcome, train_treatment, train_covariates, test_outcome, test_treatment, test_covariates, knob='dragonnet',
              output_base_dir='', binary_flag = False):
    output_dir = os.path.join(output_base_dir, knob)
    if knob == 'dragonnet':
        ate = run_data(train_outcome, train_treatment, train_covariates, test_outcome, test_treatment, test_covariates, dragon='dragonnet', binary_flag = binary_flag)

    if knob == 'tarnet':
        ate = run_data(train_outcome, train_treatment, train_covariates, test_outcome, test_treatment, test_covariates, dragon='tarnet', binary_flag = binary_flag)

    return ate