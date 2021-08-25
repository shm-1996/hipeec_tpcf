from header import *
from Galaxy import Galaxy
from Plot_Class import myPlot


def pipeline_galaxy(galaxy_name,method='masked',
    outdir=None,overwrite=False,plot=False,fit=True,age_cut=1.e7):
    """
    Perform all required operations for one galaxy.
    Parameters
        ----------
        galaxy_name : string
            Name of galaxy
        method: string
            Method to prepare random catalog.
        outdir: string
            Output directory to store summary/plots
        overwrite: Boolean
            Flag to overwrite summary file if it already exists  
        plot : Boolean 
            Flag to plot TPCF individually + other diagnostic plots 
        fit : Boolean
            Flag to perform MCMC fits   
        age_cut : Float
            The age cut to separate young and old clusters  
        Returns
        -------
        galaxy_class : Class Galaxy
            Instance of class galaxy is returned
    """
    import itertools
    print("Performing TPCF pipeline for galaxy {}.".format(galaxy_name))
    
    #Initialise galaxy class
    galaxy_class = Galaxy(galaxy_name,verbose=True)

    #Change output directory if passed by user
    if(outdir is not None):
        galaxy_class.outdir = outdir
        
    if(method not in ['masked','uniform','masked_radial','exponential']):
        raise myError("Method not recognised.")
    
    #Make sure path exists else create
    if(not os.path.exists(galaxy_class.outdir)):
        print("Path for output {} do not exist. Creating now"
            .format(galaxy_class.outdir))
        os.makedirs(galaxy_class.outdir)

     #First check if pickle file of class already exists here
    if(overwrite is False):
        if(os.path.isfile(galaxy_class.outdir+
            '/{}_summary.pkl'.format(galaxy_class.name))):
            print("Summary file exists and overwrite=False. Reading existing file.")
            galaxy_class = loadObj(galaxy_class.outdir+
            '/{}_summary'.format(galaxy_class.name))
            computeTPCF = False
        else :
            print("Summary file {}/{}_summary.pkl not present in output directory . \
Computing TPCF now..".format(galaxy_class.outdir,galaxy_class.name))  
            computeTPCF = True  
    else :
        print("Overwrite flag is provided. Computing and overwriting TPCF even if it exists.")
        computeTPCF = True

    

    #Pipeline for each galaxy starts now - 

    #TODO: Implement this

    ###################################################################################
    # 1. Compute TPCF for -
    if(method == 'exponential'):
        #Extract old TPCF r_c fitted
        sampler = loadObj(masked_directory+'MCMC_sampler')
        samples = sampler.flatchain
        PF_median_fit = np.percentile(samples,50,axis=0)
        r_c = PF_median_fit[2]
        galaxy_class.rc = r_c

    # a. all clusters
    if(computeTPCF):
        print("TPCF of all clusters.....")
        galaxy_class.Compute_TPCF(random_method=method,age=None)
        # # b. young clusters (T<10 Myr)
        print("TPCF for clusters with age <= {}".format(age_cut))
        galaxy_class.Compute_TPCF(random_method=method,age='young',age_cut=age_cut)
        # # c. old clusters (T>10 Myr)
        print("TPCF for clusters with age > {}".format(age_cut))
        galaxy_class.Compute_TPCF(random_method=method,age='old',age_cut=age_cut)
        saveObj(galaxy_class,galaxy_class.outdir+'/{}_summary'.format(galaxy_class.name))

    #Plot TPCF if required
    if(plot):
        pl = myPlot(galaxy_class)
        pl.plot_TPCF(save=True,function=None,age='both',omega1=True,
            filename=pl.galaxy.outdir+'/{}_{}Myr.pdf'.format(galaxy_name,age_cut/1.e6))
        pl.plot_TPCF(save=True,function=None,age=None,omega1=True,
            filename=pl.galaxy.outdir+'/{}_All.pdf'.format(galaxy_name))

    if(fit):

        ###################################################################################
        # 3. Fit TPCF with MCMC - both young and old all 3 functional forms
        for func,age in itertools.product(['singlepl','piecewise','singletrunc'],
            ['young','old']):
            galaxy_class.fit_TPCF(method='mcmc',function=func,age=age,plot=plot)
        saveObj(galaxy_class,galaxy_class.outdir+'/{}_summary'.format(galaxy_class.name))

        ###################################################################################
        # 4. Compare AIC for 3 models
        galaxy_class.compare_AIC(age='young')
        galaxy_class.compare_AIC(age='old')

        ###################################################################################
        # 5. Compute Inferred Physical properties: lcorr, D2 and ExpDiskrc
        #galaxy_class.getPhysicalProps()

        print("Saving class object of {} as pickle file.".format(galaxy_class.name))
        saveObj(galaxy_class,galaxy_class.outdir+'/{}_summary'.format(galaxy_class.name))
        print("\n\n")

    print("TPCF Pipeline Done..............")
    
    

    return galaxy_class

if __name__ == "__main__":

    #Parsing Arguments
    ############################################################################
    ap = argparse.ArgumentParser(description=
        'Command Line Inputs for tpcf-starclusters. All inputs optional. ')
    ap.add_argument('-method',action='store',type=str,default='masked',
        help='Method to prepare the random catalog: "Uniform","Masked"' +
        '" Masked_radial (default)" ')
    ap.add_argument('-galaxy',action='store',type=str,default=None,
        help = 'Galaxy for which tpcf to be computed. By default done for all.')
    ap.add_argument('-outdir',action='store',type=str,default=None,
        help = 'Alternate output directory for plots and files.')
    ap.add_argument('-fit',action='store_true',
        help='Flag to attempt fits with an MCMC." ')
    ap.add_argument('-plot',action='store_true',
        help='Flag to plot the TPCF." ')
    ap.add_argument('-overwrite',action='store_true',
        help='Overwrite saved pickle file." ')
    ap.add_argument('-age_cut',action='store',type=float,default=1.e7,
        help='The age cut to separate young and old clusters')
    args = vars(ap.parse_args())

    method = args['method'].lower()
    if(method not in ['uniform','masked','masked_radial','exponential']):
        raise ValueError("This method does not exist. Allowed values are "+
            "'Uniform', 'Masked', and 'Masked_Radial'.")

    galaxy_name = args['galaxy'].upper()
    if(args['outdir'] is not None):
        output_directory = os.path.abspath(args['outdir']+'/')
        if(not os.path.exists(output_directory)):
            os.makedirs(output_directory)

    else:
        output_directory = None
    #Arguments parsed
    ############################################################################
    #Only for one galaxy
    if(galaxy_name is not None):
         #for all galaxies
        if(galaxy_name == 'ALL'):
            for gname in list_of_galaxies:
                pipeline_galaxy(gname,method=method,outdir=output_directory,
                    overwrite=args['overwrite'],plot=args['plot'],fit=args['fit'],age_cut=args['age_cut'])
        elif(galaxy_name in list_of_galaxies):

            print("Running tpcf-starclusters for {}.".format(galaxy_name))
            pipeline_galaxy(galaxy_name,method=method,outdir=output_directory,
                    overwrite=args['overwrite'],plot=args['plot'],fit=args['fit'],age_cut=args['age_cut'])

        else:
            raise myError("The provided galaxy {} is not in the list of galaxies".format(galaxy_name)+
                " for which cluster catalogs are available with HIPEEC.")

   

    else :
        raise myError("Please provide galaxy name or ALL for all galaxies.")







    
        
        



