from header import *
from Plot_Class import myPlot, bbox
#from CreateTables import compare_AIC, Get_Cutoff_Scale, compare_AICc
import matplotlib.lines as lines
from regions import read_ds9
from astropy.coordinates import SkyCoord
import regions
from TPCF import linear_function,linear_truncation,\
        smooth_function,onepowerlaw_function,piecewise_truncation,\
        bootstrap_two_point,bootstrap_two_point_angular

#Axes limits in parsec
global_axes_limits = [8,1.e4]
#from PaperPlots import get_separations, deproject_region_centre
from matplotlib.patches import Patch



def Fig1(save=True,outdir='../plots/',method='masked',indir=None):

    #Create figure and axs instance
    #fig,axs = plt.subplots(nrows=4,ncols=3,figsize=(12,16))
    fig = plt.figure(figsize=(14,8),constrained_layout=True)

    if(indir == None):
        indir = os.path.abspath('../results/')+'/'
        method_dir = ''
    else :
        method_dir = None

    #Directories
    indir = os.path.abspath(indir)+'/'
    outir = os.path.abspath(outdir)+'/'

    i,j = 0,0
    #Loop through the galaxies
    for galaxy_name in list_of_galaxies:
        galaxy_dir = indir

        galaxy_class = loadObj(galaxy_dir+galaxy_name+'_summary')
        galaxy_class.get_ra_dec()
        pl = myPlot(galaxy_class)
        
        hdu = fits.open(pl.galaxy.fits_file)[1]
        wcs_galaxy = WCS(hdu.header)
        xpix,ypix = wcs_galaxy.all_world2pix(pl.galaxy.ra_raw,pl.galaxy.dec_raw,0)
        ages = galaxy_class.get_cluster_ages()
        image = hdu.data

        axs = plt.subplot2grid((2,3),(i,j),projection=wcs_galaxy,fig=fig)       

        cmap = cmr.fusion_r

       
        im1 = axs.scatter(xpix,ypix,c=np.log10(ages),
            alpha=0.8,cmap=cmap,lw=0.5,s=8.0,edgecolors='none',
            vmin=6,vmax=np.log10(np.max(ages)))
        cbar = fig.colorbar(im1,ax = axs,use_gridspec=False,
                    orientation='vertical',pad=0.00,aspect=30)
        imin,imax = axs.get_xlim()
        jmin,jmax = axs.get_ylim()

        image[np.where(image <= 1.e-2)] = np.nan
        with np.errstate(divide='ignore', invalid='ignore'):
            im = axs.imshow(np.log10(image),vmin=-2.0,alpha=0.3)

        #imin,imax,jmin,jmax = bbox(~np.isnan(image))
        xextent,yextent = (imax-imin),(jmax-jmin)
        imin,imax = imin - xextent*0.1,imax + xextent*0.1
        jmin,jmax = jmin - yextent*0.1,jmax + yextent*0.1
        axs.set_xlim(imin,imax)
        axs.set_ylim(jmin,jmax)

        ra = axs.coords['RA']
        dec = axs.coords['DEC']
        ra.set_ticklabel(exclude_overlapping=True)
        dec.set_ticklabel(exclude_overlapping=True)
        ra.set_auto_axislabel(False)
        dec.set_auto_axislabel(False)

        # #Draw 50 arcsec scale bar
        
        #No of pixels in axes
        total_pixels = imax-imin
        xmin,ymin = wcs_galaxy.all_pix2world(imin,jmin,0)
        xmax,ymax = wcs_galaxy.all_pix2world(imax,jmax,0)

        #This is in degrees
        length_per_pixel = (xmax-xmin)/(total_pixels)
        #Convert to arcsec
        length_per_pixel = length_per_pixel*3600.
        #Convert to parsec 
        length_per_pixel = pl.sep_to_pc(length_per_pixel)

        length = pl.sep_to_pc(10)
        no_pixels = np.abs(length/length_per_pixel)
        no_pixels = no_pixels/total_pixels

       
        scale = lines.Line2D([0.8,0.8+no_pixels],[0.1],
                                     lw=1,color='black',
                                    transform=axs.transAxes)
        axs.add_line(scale)
        axs.annotate(r'$10^{\prime \prime} = %d \, \mathrm{pc}$'%length,(0.7,0.15),
            xycoords='axes fraction',
                            fontsize=12)

        if(j==2):
            cbar.ax.set_ylabel(r"$\log_{10} \, \mathrm{Age} \, (\mathrm{yr})$",
                rotation=90,labelpad=5,fontsize=20)
        if(j==0):
            axs.set_ylabel(r"$\mathrm{Declination \; (J2000)}$")
        else:
            dec.set_axislabel('')

        if(i == 1):
            axs.set_xlabel(r"$\mathrm{Right \; Ascension \; (J2000)}$")
        else:
            ra.set_axislabel('')            

        axs.text(0.07,0.9,r'$\mathrm{NGC}$'+' '+r'${}$'.format(galaxy_name.split('_')[1]),
            transform=axs.transAxes)

        #Get position of subplot
        j +=1
        if(j==3):
            j = 0
            i +=1


    if(save):    
        filename = outdir+'Fig1.pdf'
        plt.savefig(filename,bbox_inches='tight')
        plt.close()
    else :
        plt.show()


    return



def Figure2(save=False,outdir='../plots/',indir=None,method='masked'):
    """
    Plot the TPCF of clusters based on age cuts.

    Parameters:
        save: 
            Flag to save the plot
        outdir: 
            Output directory in which to store plot. Default is results directory.
        indir :
            Input directory from which to read the results. Default is results directory.
        method : 
            Method for which TPCF's have been computed.   
    Returns:
        None


    """

    from matplotlib.lines import Line2D
    from regions import read_ds9
    print("Plotting TPCFs with Age cuts.")
    Merger_Stages = [5,5,6,3,5,5]

    global_axes_limits = [20,2.e4]

    #Create figure and axs instance
    fig,axs = plt.subplots(nrows=2,ncols=3,figsize=(14,8))
    i,j = 0,0

    indir = os.path.abspath('../results/')+'/' 
    outir = os.path.abspath(outdir)+'/'       

    #Directories
    indir = os.path.abspath(indir)+'/'
    for galaxy_name in list_of_galaxies:
        galaxy_class = loadObj(indir+galaxy_name+'_summary')
        plot_class = myPlot(galaxy_class)
        if(i==0):
            plot_class.plot_TPCF(function=None,age='both',axs=axs[i,j],sec_axis=True)
        else:
            plot_class.plot_TPCF(function=None,age='both',axs=axs[i,j],sec_axis=False)
        
        distance = galaxy_class.distance*const.Parsec*1.e6
        axs_limits = [0.0,0.0]
        axs_limits[0] =  global_axes_limits[0]*const.Parsec/distance*u.radian.to(u.arcsec)
        axs_limits[1] =  global_axes_limits[1]*const.Parsec/distance*u.radian.to(u.arcsec)
        axs[i,j].set_xlim(axs_limits[0],axs_limits[1])

        axs[i,j].set_ylim(0.1,1.2e2)
        
        #Figure out edge effect region
        region = read_ds9(galaxy_class.region_file)
        hdu = fits.open(galaxy_class.fits_file)[1]
        wcs = WCS(hdu.header)
        name = galaxy_class.name.split('_')[0] + ' ' + galaxy_class.name.split('_')[1]
        ra_dec = SkyCoord.from_name(name)
        xpix_c,ypix_c = wcs.all_world2pix(ra_dec.ra.value,ra_dec.dec.value,0)

        for k in range(0,np.size(region)):
            region_dep = deproject_region_centre(region,k,xpix_c,ypix_c,galaxy_class)
            if(region_dep is not None):
                sides = region_dep.to_sky(wcs).vertices
                sizes = get_separations(sides,plot_class)

        ########################################
        #Get probable boundary scale
        #TODO: Can improve this definition of scale
        boundary_scale = np.max(sizes)/5.0
        boundary_scale = plot_class.pc_to_sep(boundary_scale)
        
        axs[i,j].axvspan(boundary_scale,axs_limits[1],alpha=0.3,color='#8D717490')
        
        if(i == 0):
            axs[i,j].set_xlabel('')
        if(j>0):
            axs[i,j].set_ylabel('')
            
        if(i+j>0):
            axs[i,j].legend().remove()
        axs[i,j].text(0.1,0.1,r'$\mathrm{NGC}$'+' '+r'${}$'.format(galaxy_name.split('_')[1]),
                transform=axs[i,j].transAxes)
        
        j +=1
        if(j==3):
            j = 0
            i +=1
    if(save):    
        filename = outdir+'Fig2.pdf'
        plt.savefig(filename,bbox_inches='tight')
        plt.close()
    else :
        plt.show()


    return


def Figure3(save=False,outdir='../plots/',indir=None,method='masked'):
    """
    Plot the TPCF of clusters based on age cuts.

    Parameters:
        save: 
            Flag to save the plot
        outdir: 
            Output directory in which to store plot. Default is results directory.
        indir :
            Input directory from which to read the results. Default is results directory.
        method : 
            Method for which TPCF's have been computed.   
    Returns:
        None


    """

    fig,axs = plt.subplots(nrows=2,ncols=3,figsize=(14,8))
    i,j = 0,0
    Merger_Stages = [5,5,6,3,5,5]

    for galaxy in list_of_galaxies:
        galaxyObj = Galaxy(galaxy)
        
        if(i == 0):
            sec_axis = True

        #T<2 Myr
        galaxyObj.Compute_TPCF(age='young',age_cut=2.e6)
        pl = myPlot(galaxyObj)
        pl.plot_TPCF(age='young',axs=axs[i,j],sec_axis=sec_axis,c='#1B98F5',lw=1.0,fmt='.-')

        #T<10 Myr
        galaxyObj.Compute_TPCF(age='young',age_cut=1.e7)
        pl = myPlot(galaxyObj)
        pl.plot_TPCF(age='young',axs=axs[i,j],sec_axis=False,c='#4AA83B',lw=1.0,fmt='.-')

        #T<50 Myr
        galaxyObj.Compute_TPCF(age='young',age_cut=5.e7)
        pl = myPlot(galaxyObj)
        pl.plot_TPCF(age='young',axs=axs[i,j],sec_axis=False,c='#F5623E',lw=1.0,fmt='.-')
        
        if(i == 0):
            axs[i,j].set_xlabel('')
        if(j>0):
            axs[i,j].set_ylabel('')
                
        if(i+j>0):
            axs[i,j].legend().remove()
        axs[i,j].text(0.1,0.1,r'$\mathrm{NGC}$'+' '+r'${}$'.format(galaxy.split('_')[1]),
               transform=axs[i,j].transAxes)

        j +=1
        if(j==3):
            j = 0
            i +=1
        
    fig.savefig('../plots/AgeCuts.pdf',bbox_inches='tight')



