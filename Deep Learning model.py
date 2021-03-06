import tensorflow as tf
import numpy as np
import matplotlib.pyplot as plt


t_train = np.linspace(0, 5, 500, endpoint=True)
t_tr = np.zeros((len(t_train),1))
for i in range(len(t_train)):
  t_tr[i] = t_train[i]


x0=1
y0=1
h1=10


tf.compat.v1.disable_eager_execution()
t = tf.compat.v1.placeholder(tf.float32, [None, 1])



W1 = tf.Variable(tf.zeros([1, h1]))
b1 = tf.Variable(tf.zeros([h1]))
x1 = tf.nn.sigmoid(tf.matmul(t, W1) + b1)
W2 = tf.Variable(tf.zeros([h1, 1]))
b2 = tf.Variable(tf.zeros([1]))
NN = tf.matmul(x1, W2) + b2
x = x0 + tf.multiply(t, NN)
dx = tf.gradients(x, t)
_W1 = tf.Variable(tf.zeros([1, h1]))
_b1 = tf.Variable(tf.zeros([h1]))
y1 = tf.nn.sigmoid(tf.matmul(t, _W1) + _b1)
_W2 = tf.Variable(tf.zeros([h1, 1]))
_b2 = tf.Variable(tf.zeros([1]))
_NN = tf.matmul(y1, _W2) + _b2
y = y0 + tf.multiply(t, _NN)
_dx = tf.gradients(y, t)






#Viet ham fx gx
f = (MCBlowAir() + MCExtAir() + MCPadAir(CO2_Air = x) - MCAirTop(CO2_Air
= x, CO2_Top = y) - MCAirCan(CO2_Air = x) - MCAirOut(CO2_Air = x))
# / cap_CO2_Air
g = (MCAirTop(CO2_Air = x, CO2_Top = y) - MCTopOut(CO2_Top = y))
# / cap_CO2_Top


loss1 = tf.reduce_mean((dx - f)**2 + (_dx - g)**2)