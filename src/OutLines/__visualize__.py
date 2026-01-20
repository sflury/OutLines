from OutLines.__funcs__ import *
import matplotlib.pyplot as plt
from matplotlib.patches import Arc,Rectangle
from matplotlib import rcParams
from cycler import cycler
# vygrboi color scheme -- github.com/sflury/vygrboi
rcParams['axes.prop_cycle'] = cycler('color', ['#64008f','#fdca40','#08772f',\
                                    '#ab0014','#008cf7','#da6600','#001e52'])
rcParams['figure.figsize']       = (5,4.5) # figure size in inches
rcParams['xtick.top']            = True    # draw ticks on the top side
rcParams['xtick.bottom']         = True    # draw ticks on the bottom side
rcParams['xtick.major.size']     = 10      # major tick size in points
rcParams['xtick.minor.size']     = 4       # minor tick size in points
rcParams['xtick.major.width']    = 1       # major tick width in points
rcParams['xtick.direction']      = 'in'    # direction'] {in, out, inout}
rcParams['xtick.minor.visible']  = True    # visibility of minor ticks on x-axis
rcParams['ytick.left']           = True    # draw ticks on the left side
rcParams['ytick.right']          = True    # draw ticks on the right side
rcParams['ytick.major.size']     = 10      # major tick size in points
rcParams['ytick.minor.size']     = 4       # minor tick size in points
rcParams['ytick.major.width']    = 1       # major tick width in points
rcParams['ytick.direction']      = 'in'    # direction'] {in, out, inout}
rcParams['ytick.minor.visible']  = True    # visibility of minor ticks on y-axis


'''
Name:
    PlotIsoContours

Purpose:
    Plot isovelocity contours for a given velocity field and velocity index.

Keyword Arguments:
    :VelocityIndex (*float*): velocity field power index, default = 1
    :VelocityField (*str*): the desired velocity field, default = \'BetaCAK\'

Attributes:
    :fig (*maplotlib.pyplot.figure*): figure instance
    :ax  (*matplotlib.pyplot.axis*): axis instance
'''
class PlotIsoContours(object):
    def __init__(self,VelocityIndex=1,VelocityField='BetaCAK'):
        fig,ax = self.isovel(beta=VelocityIndex,VelocityField=VelocityField)
        self.fig = fig
        self.ax  = ax
    def isovel(self,beta=1,VelocityField='BetaCAK'):
        rscl = 10
        figure,axis = plt.subplots(1,1)
        axis.annotate('',xy=(1.02,0),xytext=(1/rscl-0.01,0),arrowprops={'arrowstyle':'->'},zorder=4)

        axis.plot([0,0],[-1,-1/rscl],'-',lw=1,color='k')
        axis.plot([0,0],[ 1/rscl, 0.2],'-',lw=1,color='k')
        axis.plot([0,0],[ 0.900,  1.0],'-',lw=1,color='k')
        axis.text(-0.05,0.26,'plane of sky',rotation=270)
        axis.text(1.07,-0.1,'to\nobserver',horizontalalignment='center')

        theta = linspace(0,2*pi,1001)
        plt.plot(cos(theta)/rscl,sin(theta)/rscl,color='k')
        for u,color in zip([-0.5,-0.3,-0.1,0.1,0.3,0.5],['C3','C5','C1','C2','C4','C0']):
            if hasattr(beta,'__len__'):
                for b,ls in zip(beta,[':','-','-.','--']):
                    uw = linspace(abs(u),v[VelocityField](rscl*0.85,b),1001)
                    rbnd = x[VelocityField](uw,b)/rscl
                    tr = arccos(abs(u)/uw)
                    xt = rbnd*cos(tr)
                    yt = rbnd*sin(tr)
                    axis.plot(sign(u)*xt+u, yt,color=color,ls=ls,lw=2)
                    axis.plot(sign(u)*xt+u,-yt,color=color,ls=ls,lw=2)
            else:
                uw = linspace(abs(u),v[VelocityField](rscl*0.9,beta),1001)
                rbnd = x[VelocityField](uw,beta)/rscl
                tr = arccos(abs(u)/uw)
                xt = rbnd*cos(tr)
                yt = rbnd*sin(tr)
                axis.plot(sign(u)*xt+u, yt,color=color,ls='-',lw=2)
                axis.plot(sign(u)*xt+u,-yt,color=color,ls='-',lw=2)
        axis.set_aspect('equal')
        axis.axis('off')
        axis.text(-1,1,VelocityField,fontsize=14,fontweight='bold')
        axis.set_xlim(-1.1,1.2)
        axis.set_ylim(-1.1,1.2)
        plt.subplots_adjust(left=0,right=1,top=1,bottom=0)
        return figure,axis

'''
Name:
    PlotGeometry

Purpose:
    Plot an edge-on visualization of the geometry of a directed outflow with
    optional cavity and optional disk.

Arguments:
    :i (*float*): inclination of the outflow with respect to the line of
        sight; if i > 6.3, interpreted as units of degrees, otherwise as radians
    :t (*float*): opening angle of the outflow with respect to the inclination;
        if t > 6.3, interpreted as units of degrees, otherwise as radians

Keyword Arguments:
    :c (*float*):  cavity angle of the outflow with respect to the inclination;
        if t > 6.3, interpreted as units of degrees, otherwise as radians
    :xd (*float*): disk radius relative to the launch radius, must be > 1

Attributes:
    :fig (*maplotlib.pyplot.figure*): figure instance
    :ax  (*matplotlib.pyplot.axis*): axis instance
'''
class PlotGeometry(object):
    def __init__(self,i,t,c=0,xd=0):
        if i > pi/2 :
            i *= pi/180
        if t > pi/2 :
            t *= pi/180
        if c > pi/2 :
            c *= pi/180
        if c > t :
            print('Opening angle must be greater than cavity angle.')
        else:
            fig,ax = self.plot_geom(i,t,c,xd)
            self.fig = fig
            self.ax  = ax
    def plot_geom(self,i,t,c,xd):
        theta = linspace(0,2*pi,1001)
        theth = linspace(0,pi,1001)
        tcone = linspace(i-t,i+t,1001)
        figure,axis = plt.subplots(1,1,figsize=(4.5,4.5),sharex=True,sharey=True)
        axis.set_aspect('equal')
        axis.axis('off')
        axis.set_xlim(-5.1,5.1)
        axis.set_ylim(-5.1,5.1)
        axis.plot([0,0],[-20,20],'--',color='k')
        axis.text(-0.55,4,'plane\nof sky',rotation=270)
        axis.annotate('',xy=(5.1,0),xytext=(-0.2,0),arrowprops={'arrowstyle':'->'},zorder=4)
        axis.text(2,-0.5,'to observer')
        # plot outflow cone
        axis.plot([5*cos(i-t+pi),5*cos(i-t)],[5*sin(i-t+pi),5*sin(i-t)],color='k',lw=2,zorder=4)
        axis.plot([5*cos(i+t+pi),5*cos(i+t)],[5*sin(i+t+pi),5*sin(i+t)],color='k',lw=2,zorder=4)

        axis.plot([0,5*cos(i)],[0,5*sin(i)],'--',color='k') # bisecting inclination angle
        axis.add_patch(Arc((0,0),3.6,3.6,angle=0,theta1=0,theta2=i*180/pi,color='C3',lw=2))
        axis.text(2.1*cos(i/2),1.9*sin(i/2),r'$i$',color='C3',horizontalalignment='center')
        axis.add_patch(Arc((0,0),5.8,5.8,angle=0,theta1=(i-t)*180/pi,theta2=i*180/pi,color='C2',lw=1.5))
        axis.add_patch(Arc((0,0),6.2,6.2,angle=0,theta1=(i-t)*180/pi,theta2=i*180/pi,color='C2',lw=1.5))
        axis.text(3.1*cos(i-t/2),4*sin(i-t/2),r'$\theta_o$',color='C2')
        # plot hollow inner part of cone
        if c > 0 :
            axis.plot([5*cos(i-c+pi),5*cos(i-c)],[5*sin(i-c+pi),5*sin(i-c)],color='k',lw=2,zorder=4)
            axis.plot([5*cos(i+c+pi),5*cos(i+c)],[5*sin(i+c+pi),5*sin(i+c)],color='k',lw=2,zorder=4)
            axis.add_patch(Arc((0,0),4.1,4.1,angle=0,theta1=i*180/pi,theta2=(i+c)*180/pi,color='C5',lw=1.5))
            axis.add_patch(Arc((0,0),4.3,4.3,angle=0,theta1=i*180/pi,theta2=(i+c)*180/pi,color='C5',lw=1.5))
            axis.add_patch(Arc((0,0),4.5,4.5,angle=0,theta1=i*180/pi,theta2=(i+c)*180/pi,color='C5',lw=1.5))
            axis.text(2.3*cos(i+c/2),2.4*sin(i+c/2),r'$\theta_c$',color='C5')
        # plot disk
        if xd > 1 :
            axis.fill_between([-10,-sin(i)*xd,sin(i)*xd],[-cos(i)*xd,-cos(i)*xd,-cos(i)*xd],[cos(i)*xd,cos(i)*xd,-cos(i)*xd],color='0.5',alpha=0.5,edgecolor='none')
            dr = 0.1
            axis.add_patch(Rectangle((-xd,-dr),2*xd,2*dr,angle=(i-pi/2)*180/pi,rotation_point=(0,0),facecolor='w',edgecolor='k',zorder=10))
            axis.text(-0.8*sin(i)*xd,0.6*cos(i)*xd,'disk',rotation=(i-pi/2)*180/pi)
        # plot source
        axis.fill_between([-10,0],[-1,-1],[1,1],color='0.5',alpha=0.5,edgecolor='none')
        axis.fill_between(cos(theth),sin(theth+pi)[::-1],sin(theth),color='w',edgecolor='k',zorder=10)
        axis.text(-0.8,-0.15,'source',zorder=11)
        plt.subplots_adjust(left=0,right=1,top=0.95,bottom=0.05)
        return figure,axis
'''
Name:
    PlotConeProjection

Purpose:
    Plot an edge-on and face-on visualization of the directed outflow cone
    projections onto the plane of the sky, including three representative
    band contours illustrating the fraction of a given de-projected velocity
    inscribed by the projected outflow.

Arguments:
    :i (*float*): inclination of the outflow with respect to the line of
        sight; if i > 6.3, interpreted as units of degrees, otherwise as radians
    :t (*float*): opening angle of the outflow with respect to the inclination;
        if t > 6.3, interpreted as units of degrees, otherwise as radians
    :c (*float*): cavity angle of the outflow with respect to the inclination;
        if t > 6.3, interpreted as units of degrees, otherwise as radians

Attributes:
    :fig (*maplotlib.pyplot.figure*): figure instance
    :ax  (*matplotlib.pyplot.axis*): axis instance
'''
class PlotConeProjection(object):
    def __init__(self,i,t,c=0):
        if i > pi/2 :
            i *= pi/180
        if t > pi/2 :
            t *= pi/180
        if c > pi/2 :
            c *= pi/180
        if c > t :
            print('Opening angle must be greater than cavity angle.')
        else:
            fig,ax = self.plot_geom(i,t,c)
            self.fig = fig
            self.ax  = ax
    @staticmethod
    def plot_cone(axis,i,t,color='C5'):
        theta = linspace(0,2*pi,1001)
        theth = linspace(0,pi,1001)
        tcone = linspace(i-t,i+t,1001)
        axis.plot(cos(theta),sin(theta),color='k',lw=2.5)
        axis.annotate('',xy=(1.02,0),xytext=(-0.02,0),arrowprops={'arrowstyle':'->'},zorder=4)
        axis.plot(cos(tcone),sin(tcone),lw=2.5,color=color,zorder=5)
        axis.plot(cos(tcone+pi),sin(tcone+pi),lw=2.5,color=color,zorder=5)
        axis.plot([cos(i-t+pi),cos(i-t)],[sin(i-t+pi),sin(i-t)],color='k',lw=1,zorder=4)
        axis.plot([cos(i+t+pi),cos(i+t)],[sin(i+t+pi),sin(i+t)],color='k',lw=1,zorder=4)
        #
        xi = cos(i)*cos(t)
        yi = sin(i)*cos(t)
        xp  = sin(t)
        yp  = 0.1*cos(i)*xp
        A  = i - 3*pi/2
        xt1 = xp*cos(theth)*cos(A)-yp*sin(theth)*sin(A)
        yt1 = xp*cos(theth)*sin(A)+yp*sin(theth)*cos(A)
        xt2 = xp*cos(theth+pi)*cos(A)-yp*sin(theth+pi)*sin(A)
        yt2 = xp*cos(theth+pi)*sin(A)+yp*sin(theth+pi)*cos(A)
        axis.plot(xt1+xi,yt1+yi,color='C5',lw=2.5,zorder=5)
        axis.plot(xt1-xi,yt1-yi,color='C5',lw=2.5,zorder=5)
        axis.plot(xt2+xi,yt2+yi,color='C5',lw=2.5,ls=':',zorder=3)
        axis.plot(xt2-xi,yt2-yi,color='C5',lw=2.5,ls=':',zorder=3)
        return axis
    @staticmethod
    def plot_band(axis,i,t,u,w,c=0,color='C0'):
        # angle of velocity projection
        rbnd = sqrt(1-(u/w)**2)
        # outer cone projected as ellipse
        f0 = sin(i)*cos(t)
        g0 = sin(t)
        h0 = g0*cos(i)
        # intersection points of band contour with outer cone projections
        a0 = (g0/h0)**2-1
        b0 = -2*f0*(g0/h0)**2
        c0 = rbnd**2+((f0/h0)**2-1)*g0**2
        q0 = (-b0+array([-1,1])*sqrt(b0**2-4*a0*c0))/(2*a0)
        # inner cone projected as ellipse
        f1 = sin(i)*cos(c)
        g1 = sin(c)
        h1 = g1*cos(i)
        # intersection points of band contour with inner cone projections
        a1 = (g1/h1)**2-1
        b1 = -2*f1*(g1/h1)**2
        c1 = rbnd**2+((f1/h1)**2-1)*g1**2
        q1 = (-b1+array([-1,1])*sqrt(b1**2-4*a1*c1))/(2*a1)
        if c == 0:
            q1 = [0,0]

        theta = linspace(-pi/2,3*pi/2,1001)
        xt = rbnd*cos(theta)
        yt = rbnd*sin(theta)
        if i == pi/2 :
            ind = where(((yt>-q0[0])&(yt<q0[0]))|(yt<-q1[0])|(yt>q1[0]))[0]
        elif i + t > pi/2 :
            ind = where(((yt>q1[0])&(yt<q1[1]))|((yt>-q0[1])&(yt<q0[0])))[0]
        elif i == 0 :
            if rbnd < h0-f0 :
                ind = array([False for j in range(len(yt))])
            else:
                ind = array([True  for j in range(len(yt))])
        else:
            ind = where(((yt>q1[0])&(yt<q1[1]))|(yt>q0[1])|(yt<q0[0]))[0]
        xt[ind] = nan
        yt[ind] = nan
        axis.plot(0.1*xt[501:]+u/w,yt[501:],zorder=4,color=color)
        axis.plot(0.1*xt[:501]+u/w,yt[:501],zorder=4,color=color,ls=':')
        xt[where(isfinite(xt))[0]] = nan
        yt[where(isfinite(yt))[0]] = nan
        xt[ind] = rbnd*cos(theta[ind])
        yt[ind] = rbnd*sin(theta[ind])

        axis.plot(0.1*xt[501:]+u/w,yt[501:],zorder=4,color=color,alpha=0.2)
        axis.plot(0.1*xt[:501]+u/w,yt[:501],zorder=4,color=color,alpha=0.2,ls=':')

        return axis
    @staticmethod
    def proj_cone(axis,i,t,color='C5'):
        theta = linspace(0,2*pi,1001)
        axis.plot(cos(theta),sin(theta),'-k',lw=2.5,zorder=0)
        if t >= pi/2 :
            axis.plot(cos(theta),sin(theta),'C5',lw=2.5)
            return axis
        # cone projected as ellipse
        f = sin(i)*cos(t)
        g = sin(t)
        h = g*cos(i)
        # intersection with unit sphere
        a = (g/h)**2-1
        b = -2*f*(g/h)**2
        xc = 1+((f/h)**2-1)*g**2
        xq = (-b+array([-1,1])*sqrt(b**2-4*a*xc))/(2*a)
        xp = sqrt(1-xq**2)
        de = arccos(xq[0])
        # flat lines at i = 90 degrees
        if i == pi/2:
            axis.plot([-g,g],[h+f,h+f],lw=2.5,color='C5')
            axis.plot([-g,g],[h-f,h-f],lw=2.5,color='C5')
            tc = linspace(pi/2-de,pi/2+de,1001)
            axis.plot(cos(tc),sin(tc),lw=2.5,color=color)
            tc = linspace(3*pi/2-de,3*pi/2+de,1001)
            axis.plot(cos(tc),sin(tc),lw=2.5,color=color)
        # if intersection with unit circle, account for lune
        # from cap projection close to transverse
        elif i + t > pi/2 :
            if isnan(xq[0]) :
                xe = g*cos(theta)
                ye = h*sin(theta)+f
                xq = 2*[ye[argmax((xe**2+ye**2))]]
                de = arccos(xq[0])
            # ellipse + lune
            xe = g*cos(theta)
            ye = h*sin(theta)+f
            xe[ye>xq[0]] = nan
            ye[ye>xq[0]] = nan
            tc = linspace(pi/2-de,pi/2+de,1001)
            axis.plot(cos(tc),sin(tc),lw=2.5,color=color)
            axis.plot(xe,ye,color='C5',lw=2.5)
            # lune
            xe = g*cos(theta)
            ye = h*sin(theta)-f
            xe[ye>-xq[1]] = nan
            ye[ye>-xq[1]] = nan
            tc = linspace(3*pi/2-de,3*pi/2+de,1001)
            axis.plot(cos(tc),sin(tc),lw=2.5,color=color)
            axis.plot(xe,ye,color='C5',lw=2.5)
        # if no intersection, just an ellipse
        else:
            xe = g*cos(theta)
            ye = h*sin(theta)+f
            axis.plot(xe,ye,color='C5',lw=2.5)
        return axis
    @staticmethod
    def proj_band(axis,i,t,u,w,c=0,color='C0'):
        # angle of velocity projection
        rbnd = sqrt(1-(u/w)**2)
        # outer cone projected as ellipse
        f0 = sin(i)*cos(t)
        g0 = sin(t)
        h0 = g0*cos(i)
        # intersection points of band contour with outer cone projections
        a0 = (g0/h0)**2-1
        b0 = -2*f0*(g0/h0)**2
        c0 = rbnd**2+((f0/h0)**2-1)*g0**2
        q0 = (-b0+array([-1,1])*sqrt(b0**2-4*a0*c0))/(2*a0)
        # inner cone projected as ellipse
        f1 = sin(i)*cos(c)
        g1 = sin(c)
        h1 = g1*cos(i)
        # intersection points of band contour with inner cone projections
        a1 = (g1/h1)**2-1
        b1 = -2*f1*(g1/h1)**2
        c1 = rbnd**2+((f1/h1)**2-1)*g1**2
        q1 = (-b1+array([-1,1])*sqrt(b1**2-4*a1*c1))/(2*a1)
        if c == 0:
            q1 = [0,0]

        theta = linspace(-pi/2,3*pi/2,1001)
        xt = rbnd*cos(theta)
        yt = rbnd*sin(theta)
        axis.plot(xt,yt,color=color,alpha=0.2)
        if i == pi/2 :
            ind = where(((yt>-q1[0])&(yt<q1[0]))|(yt<-q0[0])|(yt>q0[0]))[0]
        elif i + t > pi/2 :
            if i + c > pi/2 : #q1[1] >= 1 :
                ind = where((yt<-q1[1])|(yt>q1[0])|((yt>-q0[1])&(yt<q0[0])))[0]
            else:
                ind = where(((yt>q1[0])&(yt<q1[1]))|((yt>-q0[1])&(yt<q0[0])))[0]
        elif q0[0] < 0 :
            ind = where(((yt>q1[0])&(yt<q1[1]))|(yt<q0[0]))[0]
        elif i == 0 :
            if rbnd < h0-f0 :
                ind = array([False for j in range(len(yt))])
            else:
                ind = array([True  for j in range(len(yt))])
        else:
            ind = where(((yt>q1[0])&(yt<q1[1]))|(yt>q0[1])|(yt<q0[0]))[0]
        xt[ind] = nan
        yt[ind] = nan
        axis.plot(xt,yt,color=color,zorder=0)

        return axis
    def plot_geom(self,i,t,c=0):
        # convert to degrees
        istr = f'{int(round(i*180/pi)):0>2d}'
        tstr = f'{int(round(t*180/pi)):0>2d}'
        # plotting
        fig,ax = plt.subplots(1,2,figsize=(7,4),sharex=True,sharey=True)
        #plot the cones
        ax[0] = self.plot_cone(ax[0],i,t,color='C5')
        ax[1] = self.proj_cone(ax[1],i,t,color='C5')
        if c > 0 :
            ax[0] = self.plot_cone(ax[0],i,c,color='k')
            ax[1] = self.proj_cone(ax[1],i,c,color='k')
        ax[0].set_title('Edge-On View')
        ax[1].set_title('Plane of Sky')
        for uj,cj in zip([0.3,0.45,0.6],['C2','C4','C0']):
            # band contours
            ax[0] = self.plot_band(ax[0],i,t,uj,0.75,color=cj,c=c)
            ax[1] = self.proj_band(ax[1],i,t,uj,0.75,color=cj,c=c)
        # label
        fig.text(0.45,0.9,r'$i\ =$'+f'{int(round(i*180/pi)): >3d}'+r'$^\circ$')
        fig.text(0.45,0.825,r'$\theta_o\ =$'+f'{int(round(t*180/pi)): >3d}'+r'$^\circ$')
        if c != 0 :
            fig.text(0.45,0.75,r'$\theta_c\ =$'+f'{int(round(c*180/pi)): >3d}'+r'$^\circ$')
        #ax[1].text(0,0,r'$\oplus$',color='k')
        for a in ax:
            a.set_aspect('equal')
            a.axis('off')
            a.set_xlim(-1.1,1.1)
            a.set_ylim(-1.1,1.1)
        plt.subplots_adjust(left=0,right=1,top=0.95,bottom=0)
        return fig,ax
