# QPSK signaling simulation
import numpy as np
import matplotlib.pyplot as plt

def main():
    # source input
    t1 = np.arange(0, 8.5, 0.5)

    plt.subplot(4, 1, 1)
    y1 = [0, 1, 0, 1, 1, 0, 1, 0, 1, 1, 0, 1, 0, 1, 0, 1, 1]                              #
    plt.plot(t1, y1, drawstyle='steps-post')
    plt.xlim(0, 8)
    plt.ylim(-0.5, 1.5)
    plt.title('Input Signal')

    # common variable of I and Q channel
    a = 1 / np.sqrt(2) # normalized amplitude
    tI = np.arange(0, 9, 1)

    # Real part, I Signal
    plt.subplot(4, 1, 2)
    yI = [-a, a, -a, a, -a, a, -a, a, a]
    plt.plot(tI, yI, drawstyle='steps-post')
    plt.xlim(0, 8)
    plt.ylim(-2, 2)
    plt.title('I Signal')

    # Imaginary part, Q Signal
    plt.subplot(4, 1, 3)
    yQ = [a, -a, -a, a, -a, a, -a, a, a]
    plt.plot(tI, yQ, drawstyle='steps-post')
    plt.xlim(0, 8)
    plt.ylim(-2, 2)
    plt.title('Q Signal')

    # QPSK Signal
    plt.subplot(4, 1, 4)
    t = np.arange(0, 9, 0.01)

    def output_waveform(I, Q, t):
        rect_wave = []
        for i in range(0, len(I)):
            t_temp = t[((i)*100) : ((i+1)*100)]
            yI_temp = yI[i]*np.ones(100)
            yQ_temp = yQ[i]*np.ones(100)
            wave_temp = yI_temp*np.cos(2*np.pi*5*t_temp) - yQ_temp*np.sin(2*np.pi*5*t_temp)
            rect_wave.append(wave_temp)
        return rect_wave

    rect_wave = output_waveform(yI, yQ, t)
    plt.plot(t, np.array(rect_wave).flatten(), 'r')
    plt.xlim(0, 8)
    plt.ylim(-2, 2)
    plt.title('QPSK Signal')

    plt.tight_layout()
    plt.show()

if __name__ == '__main__':
    main()
