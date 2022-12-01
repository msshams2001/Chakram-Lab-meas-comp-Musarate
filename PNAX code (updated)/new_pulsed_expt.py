__author__ = 'Alex Ma/Nate Earnest'

from slab import *
from numpy import *
import os
from tqdm import tqdm
from slab.dataanalysis import get_next_filename

im = InstrumentManager()
def polar2mag(xs, ys):
    return 20*np.log10(np.sqrt(xs ** 2 + ys ** 2)), np.arctan2(ys, xs) * 180/np.pi

class PnaxPulseExperiment(Experiment):
    '''pnax pulsed experiments'''

    def __init__(self, config_file = '..\\config_pnax.json', **kwargs):

        self.__dict__.update(kwargs)


        super(PnaxPulseExperiment, self).__init__(path=os.getcwd(), prefix='', config_file=config_file, liveplot_enabled=False, **kwargs)   
            
        print('Writing to datafile:', self.fname)

        # this formats the config file into a single line....
        # self.save_config()
        self.pnax = im['PNAX2']
        # previous line already loads the devices in 'aliases'

        print('PNAX Connected:')
        print(self.pnax.get_id())
        # print(self.dcflux.get_id())
        # print(self.drive2.get_id())

        # if self.cfg['flux']['set_flux']:
        #     self.dcflux.set_mode('current')
        #     # self.dcflux.set_range(self.cfg['flux']['range'])
        #     # self.dcflux.set_output(True)
        #     print("Dcflux range set to %.4f mA" % (1000*self.cfg['flux']['range']))
        #     print("Ramping dcflux to %.4f mA" % (1000*self.cfg['flux']['current']))
        #     self.dcflux.set_current(self.cfg['flux']['current'])

        # initialize all common pulse settings
        self.initialize_pnax()
        print('Readout power is %.2f'%self.cfg['readout']['power'])

        try:
            self.LSG = im['LB1']
            self.LSG_state=True
            print('LSG connected')
        except:
            print('Labbrick not connected, cannot do EF pulse probe')
            self.LSG_state=False

    def focus_pts(self, f0, span, steps, spacing):
            steps=int(steps/2)
            dx = 1.0 / (steps - 1)
            us = [f0 + (i * dx) ** spacing * span / 2 for i in range(steps + 1)][1::]
            ls = [f0 - (i * dx) ** spacing * span / 2 for i in range(steps)]
            ls.reverse()
            return np.array(ls + us[0:-1])

    def initialize_pnax(self):

        print("\nConfiguring the PNAX for pulse measurements")
        self.pnax.set_timeout(10E3)
        self.pnax.clear_traces()
        self.pnax.setup_measurement("S21")
        # right now ifbw cover the whole readout period, so average is done using set_sweep_points
        avgs = 1
        self.pnax.setup_take(averages_state=True)
        self.pnax.set_averages_and_group_count(avgs, True)
        # setting average - will be reset in each exp function
        self.pnax.set_sweep_points(100)

        self.pnax.write('SENSE:FOM:STATE 1')
        # todo: maybe set al ranges to uncoupled and setting freq manually??
        self.pnax.write("sense:fom:range2:coupled 0") # source 1
        self.pnax.write("sense:fom:range3:coupled 0") # receivers
        self.pnax.write("sense:fom:range4:coupled 1") # source 2

        # since range 4 is coupled, setting range1 freq effectively sets range 4 as well
        self.pnax.write('SENSE:FOM:RANGE1:FREQUENCY:START %f' % self.cfg['qubit']['freq'])
        self.pnax.write('SENSE:FOM:RANGE1:FREQUENCY:STOP %f' % self.cfg['qubit']['freq'])
        print('readout freq is: %.2f'%(self.cfg['readout']['freq']))

        self.pnax.write('SENSE:FOM:RANGE2:FREQUENCY:START %f' % (self.cfg['readout']['freq']))
        self.pnax.write('SENSE:FOM:RANGE2:FREQUENCY:STOP %f' % self.cfg['readout']['freq'])
        self.pnax.write('SENSE:FOM:RANGE3:FREQUENCY:START %f' % self.cfg['readout']['freq'])
        self.pnax.write('SENSE:FOM:RANGE3:FREQUENCY:STOP %f' % self.cfg['readout']['freq'])

        # setting the port powers and decoupling the powers
        self.pnax.write("SOUR:POW:COUP OFF")
        self.pnax.write(":SOURCE:POWER1 %f" % (self.cfg['readout']['power']))
        self.pnax.write(":SOURCE1:POWER1:MODE ON")
        self.pnax.write(":SOURCE:POWER3 %f" % (self.cfg['drive']['power']))
        self.pnax.write(":SOURCE1:POWER3:MODE ON")

        #Turning the leveling of the ports to 'open-loop'. This is required when trying to pulsed modulated measurements otherwise the ALC will try to level the source with the detected power level with pulse on and off, causing a source unleveled error
        self.pnax.write("source1:power1:alc:mode openloop")
        self.pnax.write("source1:power3:alc:mode openloop")
        #Setting up the proper trigger mode. Want to trigger on the point

        # ext sour is not necessary, pulse0 ensure int. sync? - alex ma
        #self.pnax.write("TRIG:SOUR EXT")
        self.pnax.write("SENS:SWE:TRIG:MODE POINT")

        #Turning off the "Reduce IF BW at Low Frequencies" because the PNA-X adjusts the BW automatically to correct for roll-off at low frequencies
        self.pnax.write("SENS:BWID:TRAC OFF")

        #turning on the pulses
        self.pnax.write("SENS:PULS0 1") # automatically sync ADC to pulse gen
        self.pnax.write("SENS:PULS1 1")
        self.pnax.write("SENS:PULS2 1")
        self.pnax.write("SENS:PULS3 1")
        self.pnax.write("SENS:PULS4 1")
        #turning off the inverting
        self.pnax.write("SENS:PULS1:INV 0")
        self.pnax.write("SENS:PULS2:INV 0")
        self.pnax.write("SENS:PULS3:INV 0")
        self.pnax.write("SENS:PULS4:INV 0")

        # setting the exp period
        self.pnax.write("SENS:PULS:PERiod %.12f" % (self.cfg['expt_trigger']['period_ns']*1e-9))

        # This is setting the ifbw to be wide enough to ensure that the measurement occurs within the readout pulse
        ifbw = 1.0/(0.99*self.cfg['readout']['width']*1e-9)
        self.pnax.set_ifbw(ifbw)
        # todo attention: only discrete values of ifbw
        print('target ifbw =', ifbw)
        print('actual ifbw =', self.pnax.get_ifbw())

        self.pnax.write("SENS:PULS1:WIDT %.12f" % (self.cfg['readout']['width']*1e-9))
        self.pnax.write("SENS:PULS2:WIDT %.12f" % (self.cfg['drive']['width']*1e-9))
        self.pnax.write("SENS:PULS3:WIDT %.12f" % (0.0)) # redundant
        self.pnax.write("SENS:PULS4:WIDT %.12f" % (0.0)) # powering Nate's box

        # delays of various pulses

        # take care of the case where readout is inside drive pulse
        # start delay of pulse 1
        self.read_delay = 0.99 * self.cfg['expt_trigger']['period_ns'] - \
                     max(self.cfg['readout']['delay'], self.cfg['readout']['width'])

        drive_end_delay_unadjust = self.read_delay - self.cfg['readout']['delay']
        # adjust for box response, ns
        self.pulse2end_delay = drive_end_delay_unadjust + self.cfg['pulse_box_delay_adjust']['pulse2']
        self.pulse3end_delay = drive_end_delay_unadjust + self.cfg['pulse_box_delay_adjust']['pulse3']
        self.pulse4end_delay = drive_end_delay_unadjust + self.cfg['pulse_box_delay_adjust']['pulse4']

        drive1_delay = self.pulse2end_delay - self.cfg['drive']['width']

        self.pnax.write("SENS:PULS0:DEL %.12f" % ( (self.read_delay + self.cfg['readout']['ADC_delay']) *1e-9) )
        self.pnax.write("SENS:PULS1:DEL %.12f" % (self.read_delay*1e-9) )
        self.pnax.write("SENS:PULS2:DEL %.12f" % (drive1_delay*1e-9) )
        self.pnax.write("SENS:PULS3:DEL %.12f" % (0) )
        self.pnax.write("SENS:PULS4:DEL %.12f" % (0) )

    def take_one_pulse_probe_point(self, scan_pt, use_pi_time=True, avgs=None):
        if avgs==None:
            self.pnax.set_sweep_points(self.cfg['pulse_probe']['averages'])
        else:
            self.pnax.set_sweep_points(avgs)

        if use_pi_time==True:
            drive_ns = self.cfg['pulse_info']['pi_time']
            self.pnax.write("SENS:PULS2:WIDT %.12f" % (drive_ns * 1e-9))
            drive1_delay = self.pulse2end_delay - drive_ns
            self.pnax.write("SENS:PULS2:DEL %.12f" % (drive1_delay * 1e-9))
            self.pnax.write(":SOURCE1:POWER3:MODE ON")
        else:
            self.pnax.write("SENS:PULS2:WIDT %.12f" % (self.cfg['drive']['width']*1e-9))
            drive1_delay = self.pulse2end_delay - self.cfg['drive']['width']
            self.pnax.write("SENS:PULS2:DEL %.12f" % (drive1_delay*1e-9) )

        self.pnax.write('SENSE:FOM:RANGE1:FREQUENCY:START %f' % scan_pt)
        self.pnax.write('SENSE:FOM:RANGE1:FREQUENCY:STOP %f' % scan_pt)
        data = self.pnax.take_one_in_mag_phase()
        mags = data[1]
        phases = data[2]
        mags_averaged = average(mags)
        phases_averaged = average(phases)
        return mags_averaged, phases_averaged

