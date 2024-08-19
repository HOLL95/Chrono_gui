import numpy as np
import matplotlib.pyplot as plt
from pandas import DataFrame
from matplotlib.widgets import Button, Slider
from scipy.interpolate import CubicSpline
class ChronoGui:
    def __init__(self, time, current, **kwargs):
        if "dumpfile" not in kwargs:
            kwargs["dumpfile"]="Saved_currents.csv"
        self.savename=kwargs["dumpfile"]
        self.colours=plt.rcParams['axes.prop_cycle'].by_key()['color'][1:]
        self.current=current
        self.time=time
        fig, self.ax = plt.subplots()
        self.ax.set_xlabel("Current")
        self.current_line,=self.ax.plot(time, current, label="Experimental current")
        self.ax.set_ylabel("Current")
        fig.set_size_inches(8, 6)
        self.fig=fig
        self.scatter_points=None
        self.mean_lines=None
        self.plot_lines=None
        self.mean_value=None
        self.saveax = fig.add_axes([0.8, 0.925, 0.1, 0.04])
        self.button = Button(self.saveax, 'Dump!', hovercolor='0.975')
        self.button.on_clicked(self.dump)
    def dump(self, event):
        if self.mean_lines is not None:
            time_chunks=np.array(self.time_chunks)
            save_dict={"Time start":time_chunks[:,0], "Time end":time_chunks[:,1], "Mean current (A)":self.mean_value}
            DataFrame(save_dict).to_csv(self.savename)
    
class Discontinuous(ChronoGui):
    def __init__(self, time, current, **kwargs):
        if "dumpfile" not in kwargs:
            super().__init__(time, current)
        else:
            super().__init__(time, current, dumpfile=kwargs["dumpfile"])
        self.diff=np.abs(np.diff(current))
        self.time_diff=np.mean(np.diff(time))
        self.mean_current_diff=np.mean(self.diff)
        self.ratio=self.diff/self.mean_current_diff
       
        self.fig.subplots_adjust(bottom=0.25)
        ratio_cutoff_ax = self.fig.add_axes([0.125, 0.15, 0.775, 0.03])
        time_between_pulse_ax= self.fig.add_axes([0.125, 0.1, 0.775, 0.03])
        pulse_cutoff_ax= self.fig.add_axes([0.125, 0.05, 0.775, 0.03])
        plt.ylabel("Current (A)")
        self.ratio_slider=Slider(
                ax=ratio_cutoff_ax,
                label='Cutoff',
                valmin=10,
                valmax=300,
                valinit=50,
            )
        self.pulse_slider=Slider(
                ax=time_between_pulse_ax,
                label='Pulse time',
                valmin=1,
                valmax=200,
                valinit=time[-1]/50,
            )
        self.exclude_slider=Slider(
                ax=pulse_cutoff_ax,
                label='Exclude no.',
                valmin=0,
                valmax=50,
                valinit=0,
                valfmt="%i"
            )

        for slider in [self.ratio_slider, self.pulse_slider, self.exclude_slider]:
            slider.on_changed(self.update)
       
        plt.show()
    def update(self, val):
        exclude=int(self.exclude_slider.val//1)
        self.exclude_slider.valtext.set_text(str(exclude))
        self.calculate(self.ratio_slider.val, self.pulse_slider.val, exclude)
        self.ax.set_xlim(self.time_min_max)
        self.fig.canvas.draw_idle()
    def calculate(self, cutoff, pulse, exclude):
        for line_collection in [self.mean_lines, self.plot_lines]:
            if line_collection is not None:
                for line in line_collection:
                    line.remove()
        if self.scatter_points is not None:
            self.scatter_points.remove()
       

        peak_loc=np.where(self.ratio>cutoff)
        approx_time_between_pulses=pulse

        peak_times=self.time[1:][peak_loc]
        time_locations=[peak_times[0]]
        global_idx=np.where(self.time>time_locations[0]-(10*self.time_diff))
        self.plot_lines=[]
        self.mean_lines=[]
       
        plot_time=self.time[global_idx]
        self.current_line.set_xdata(plot_time)
        self.current_line.set_ydata(self.current[global_idx])
        self.time_min_max=[plot_time[0]-(20*self.time_diff), plot_time[-1]+(20*self.time_diff)]
        for i in range(0, len(peak_times)):
            current_time=time_locations[-1]
            if peak_times[i]<(current_time+approx_time_between_pulses):
                continue
            else:
                time_locations.append(peak_times[i])
        time_locations=np.array(time_locations)
        self.time_chunks=np.zeros((len(time_locations)+1, 2))
        self.mean_value=np.zeros(len(time_locations)+1)
        inv_index=np.where(self.time<(time_locations[0]-(10*self.time_diff)))
        inv_time=self.time[inv_index]
        inv_current=self.current[inv_index]
        self.time_chunks[0,:]=[inv_time[0], inv_time[-1]]
        self.mean_value[0]=np.mean(inv_current)
        time=self.time

        current=self.current
        peaks=[current[np.where(time==x)][0] for x in time_locations]

        self.scatter_points=self.ax.scatter(time_locations, peaks,color="red", label="Shift points")
       
        exclude_timestep=exclude
        plot_colours=self.colours*(int(len(time_locations)//len(self.colours))+1)
        for i in range(1, len(time_locations)+1):
            
            if i==len(time_locations):
                second_idx=time[-1]
            else:
                second_idx=time_locations[i]

           
            data_chunk=np.where((time>(time_locations[i-1]+(exclude_timestep*self.time_diff))) & (time<second_idx-(exclude_timestep*self.time_diff)))
            
            mean=np.mean(current[data_chunk])
            plot_current=current[data_chunk]
            plot_time=time[data_chunk]
            if i==0:
                mean_line,=self.ax.plot(plot_time, np.ones(len(plot_current))*mean, color="black", linestyle="--", label="Extracted means")
            else:
                mean_line,=self.ax.plot(plot_time, np.ones(len(plot_current))*mean, color="black", linestyle="--")
            extracted_line,=self.ax.plot(plot_time, plot_current, color=plot_colours[i+1])
            self.mean_lines.append(mean_line)
            self.plot_lines.append(extracted_line)
            self.time_chunks[i,:]=[plot_time[0], plot_time[-1]]
            self.mean_value[i]=mean
            self.ax.legend()

class Smooth(ChronoGui):
    def __init__(self, time, current, **kwargs):
        if "dumpfile" not in kwargs:
            super().__init__(time, current)
        else:
            super().__init__(time, current, dumpfile=kwargs["dumpfile"])
        self.diff=np.abs(np.diff(current))
        self.time_diff=np.mean(np.diff(time))
        self.mean_current_diff=np.mean(self.diff)
        self.ratio=self.diff/self.mean_current_diff
       
        self.fig.subplots_adjust(bottom=0.3)
        pulse_interval = self.fig.add_axes([0.125, 0.25, 0.775, 0.03])
        pulse_time= self.fig.add_axes([0.125, 0.2, 0.775, 0.03])
        first_injection= self.fig.add_axes([0.125, 0.15, 0.775, 0.03])
        offset= self.fig.add_axes([0.125, 0.1, 0.775, 0.03])
        num_spline_points=self.fig.add_axes([0.125, 0.05, 0.775, 0.03])

        plt.ylabel("Current (A)")
        self.interval_slider=Slider(
                ax=pulse_interval,
                label='Interval time',
                valmin=10,
                valmax=1000,
                valinit=300,
            )
        self.pulse_slider=Slider(
                ax=pulse_time,
                label='Pulse time',
                valmin=1,
                valmax=300,
                valinit=100,
            )
        self.first=Slider(
                ax=first_injection,
                label='First injection',
                valmin=0,
                valmax=500,
                valinit=50,
            )
        self.offset=Slider(
                ax=offset,
                label='Offset',
                valmin=0,
                valmax=200,
                valinit=0,
            )
        self.num_points=Slider(
                ax=num_spline_points,
                label='Spline points',
                valmin=10,
                valmax=300,
                valinit=50,
                valfmt="%i"
            )

        

        for slider in [self.interval_slider, self.pulse_slider, self.first, self.offset, self.num_points]:
            slider.on_changed(self.update)
       
        plt.show()
    def update(self, val):
        num_points=int(self.num_points.val//1)
        self.num_points.valtext.set_text(str(num_points))
        self.calculate(self.interval_slider.val, self.pulse_slider.val, self.first.val, self.offset.val,num_points)
        self.fig.canvas.draw_idle()
    def calculate(self, pulse_interval, pulse_time, first_injection, offset, num_spline_points):
        for line_collection in [self.mean_lines, self.plot_lines]:
            if line_collection is not None:
                for line in line_collection:
                    line.remove()
        if self.scatter_points is not None:
            self.scatter_points.remove()
        spline_step=int(len(self.time)//num_spline_points)
        reduced_time=self.time[::spline_step]
        reduced_current=self.current[::spline_step]
        spline_interpolant=CubicSpline(reduced_time, reduced_current)
        second_derivative=spline_interpolant(reduced_time,2)
        min_values=sorted(enumerate(second_derivative), key=lambda x:x[1])
        
        time_array=[reduced_time[min_values[0][0]]]
        scatter_array=[]
        self.plot_lines=[]
        self.mean_lines=[]
        for i in range(1, len(min_values)):
            time_pos=reduced_time[min_values[i][0]]
            window=[time_pos-pulse_interval, time_pos+pulse_interval]
            new_peak=True
            if time_pos<first_injection:
                continue
            for j in range(0, len(time_array)):
                if time_array[j]>window[0] and time_array[j]<window[1]:
                    new_peak=False
                    break
            if new_peak==True:
                time_array.append(time_pos)
                
        self.time_chunks=np.zeros((len(time_array), 2))
        self.mean_value=np.zeros(len(time_array))
        plot_colours=self.colours*(int(len(time_array)//len(self.colours))+1)
        scatter_values=np.zeros((len(time_array), 2))
        time=self.time
        current=self.current
        for i in range(0, len(time_array)):
            second_idx=time_array[i]-offset
            first_idx=second_idx-pulse_time
            time_idx=np.where((time>(first_idx))&(time<=second_idx))
            plot_time=time[time_idx]
            plot_current=current[time_idx]
            line1,=self.ax.plot(plot_time, plot_current, color=plot_colours[i])
            mean_value=np.mean(plot_current)
            if i==0:            
                line2, =self.ax.plot(plot_time, np.ones(len(plot_time))*mean_value, color="black", linestyle="--", label="Mean")
            else:
                line2, =self.ax.plot(plot_time, np.ones(len(plot_time))*mean_value, color="black", linestyle="--")
            self.plot_lines.append(line1)
            self.mean_lines.append(line2)
            self.time_chunks[i,:]=[plot_time[0], plot_time[1]]
            self.mean_value[i]=mean_value
            scatter_values[i,:]=[plot_time[-1], plot_current[-1]]
        self.scatter_points=self.ax.scatter(scatter_values[:,0], scatter_values[:,1],color="red", label="Shift points")
        self.ax.legend()



