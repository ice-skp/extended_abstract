steps=5;

clustername = oshostname;
% {{{
cluster=generic('name',oshostname(),'np',4);
clear clustername
% }}}
org=organizer('repository',['./Models'],'prefix',['iceshelf_'],'steps',steps, 'color', '34;47;2'); clear steps;
setenv('CA_DEBUG_TRANSACTIONS', '2');

%step 1: create mesh

if perform(org,'Mesh')% {{{
	md=model;
	md=bamg(model,'domain','Square.exp','hmax',250);

	savemodel(org,md);
end %}}}

%step 2: set parameters 

if perform(org,'Param_decoup')% {{{

	md=loadmodel(org,'Mesh');
	md=setmask(md,'all','');

	%Geometry
	Hmin=200;
	Hmax=1400;
	Lx=max(md.mesh.x)-min(md.mesh.x);
    Ly=max(md.mesh.y)-min(md.mesh.y);
	ymin=min(md.mesh.y);
	ymax=max(md.mesh.y);
	xmin=min(md.mesh.x);
	xmax=max(md.mesh.x);

	md.geometry.thickness=Hmax+(Hmin-Hmax)*(md.mesh.y-ymin)/(ymax-ymin);
	md.geometry.base=-md.materials.rho_ice/md.materials.rho_water*md.geometry.thickness; %set as constants
	md.geometry.surface=md.geometry.base+md.geometry.thickness;
	md.geometry.bed=-900*ones(md.mesh.numberofvertices,1);

	%Materials
	md.initialization.temperature=(273-15)*ones(md.mesh.numberofvertices,1);
	md.materials.rheology_B=1.68e8*ones(md.mesh.numberofelements,1); % set as constant in Serg.
	md.materials.rheology_law = 'BuddJacka';
	md.materials.rheology_n=3*ones(md.mesh.numberofelements,1);

	%Surface mass balance and basal melting
	md.smb.mass_balance=zeros(md.mesh.numberofvertices,1);
	md.basalforcings.groundedice_melting_rate=zeros(md.mesh.numberofvertices,1); 
    md.basalforcings.floatingice_melting_rate=10*ones(md.mesh.numberofvertices,1).*(1-(md.mesh.y/Ly).^(1/3));

	%Friction
	md.friction.coefficient=20*ones(md.mesh.numberofvertices,1);
	md.friction.coefficient(find(md.mask.ocean_levelset<0.))=0.;
	md.friction.p=ones(md.mesh.numberofelements,1);
	md.friction.q=ones(md.mesh.numberofelements,1);

	%Numerical parameters
	md.masstransport.stabilization=1;

	%Deal with boundary conditions:
	md=SetIceShelfBC(md,'SquareFront.exp');

	%set velocity on inflow boundary
	pos=find(md.mesh.y==0);
	u0=0; uc=1000;
    Lx_half=(max(md.mesh.x)-min(md.mesh.x))/2;
	u=u0+(uc-u0)*(1-((md.mesh.x-Lx_half)./Lx_half).^4);
	md.stressbalance.spcvy(pos)=u(pos); 
	savemodel(org,md);
end
%}}}

%step 3: run SSA to steady state (140 yrs)

if perform(org,'SSATransient_decoup')% {{{

	md=loadmodel(org,'Param_decoup');
    %md=extrude(md,10,1.0);
	md=setflowequation(md,'SSA','all');
	md.miscellaneous.name='SSATransient_decoup';

	md.initialization.vx=zeros(md.mesh.numberofvertices,1);
	md.initialization.vy=zeros(md.mesh.numberofvertices,1);
	md.initialization.vz=zeros(md.mesh.numberofvertices,1);

	md.initialization.temperature = 273-15; % set surface temp to -15 everywhere 
	md.transient.isthermal=0;
	md.transient.issmb=0;

	md.timestepping.start_time=0;
	md.timestepping.time_step=1/12;

	md.timestepping.final_time=140; 
	%md.settings.output_frequency=12;
	%md.settings.checkpoint_frequency=12;
	md.verbose.solution = 1;
	md.stressbalance.requested_outputs={'default',...
		'StrainRatexx','StrainRateyy','StrainRatezz',...
		'StrainRatexy','StrainRateyz','StrainRatexz','StrainRateeffective',...
        'DeviatoricStresseffective','DeviatoricStressxx','DeviatoricStressxy','DeviatoricStressyy'};

	%GO solve
	md.cluster=cluster;
	md=solve(md,'tr'); 
    savemodel(org,md);

    md.results.TransientSolutionSmall=md.results.TransientSolution(end-1:end);
    md.results.TransientSolution=[];
    md.results.TransientSolution=md.results.TransientSolutionSmall;
    save('./Models/iceshelf_SSATransient_decoup_small.mat','md');

end%}}}

%Step 4: run picop for 70 yrs
% {{{ Picop_transient: 
if perform(org,'Picop_transient_decoup')

	md=loadmodel(org,'SSATransient_decoup_small');

    md=transientrestart(md,2); 
	md.timestepping.start_time=2/12;

	timesteps=1/12;
	md.inversion.iscontrol=0;
	md.timestepping.final_time=50;
	md.timestepping.time_step=timesteps;
	%md.settings.output_frequency=timesteps*12;
	%md.settings.recording_frequency=timesteps*12;

	md.transient.amr_frequency=0;
	md.transient.isgroundingline=0;
	md.transient.isthermal=0;
	md.transient.ismovingfront=0;
	md.transient.requested_outputs={'default',...
		'StrainRateeffective','StrainRatexx','StrainRatexy','StrainRateyy',...
		'StrainRatexz','StrainRateyz','StrainRatezz','BasalforcingsFloatingiceMeltingRate','BasalforcingsGroundediceMeltingRate'...
        'DeviatoricStresseffective','DeviatoricStressxx','DeviatoricStressxy','DeviatoricStressyy'};

	%Parameterized melting rate
	md.basalforcings=basalforcingspico(md.basalforcings);
	md.basalforcings.basin_id=ones(md.mesh.numberofelements,1);
	md.basalforcings.num_basins=1;
	md.basalforcings.maxboxcount=5;
	md.basalforcings.isplume=1; %1 for PICOP (PICO + Plume), 0 for PICO

	md.verbose=verbose('solution',true,'convergence',true);
	md.cluster=cluster;
	Smean=34.73;
	Tmean=272.53;

	md.basalforcings.farocean_salinity=[Smean Smean;timesteps md.timestepping.final_time]; %Time series
	md.basalforcings.farocean_temperature=[Tmean Tmean;timesteps md.timestepping.final_time];   %Time series
	md.miscellaneous.name='picop_transient_pert_decoup';

	md=solve(md,'tr');
	savemodel(org,md)
end
% }}}

%Step 5: run picop with channels
% {{{ 
if perform(org,'Picop_transient_pert_decoup')

	md=loadmodel(org,'SSATransient_decoup_small'); %load last 2 steps of steady state

    ind=find(md.mesh.y==0);
	Hmax=md.results.TransientSolution(end).Thickness(ind);

    x=0:10:50000;
	delta_H =200+200*((cos(x*pi*30/Lx))); %perturbation, translation of +200 to ensure thickness does not exceed 1400
	
	xi=md.mesh.x(ind); % finding grounding line
	interp_H = interp1(x,delta_H,xi); %interpolating perturbation onto mesh at gl
	%plot(x,delta_H,'o',xi,interp_H);

	%md.geometry.thickness=Hmax+(Hmin-Hmax)*(md.mesh.y-ymin)/(ymax-ymin);
	md.geometry.thickness(ind)=Hmax-interp_H; %needed to change size of 1400 to match that of H_pert
    md.geometry.base=-md.materials.rho_ice/md.materials.rho_water*md.geometry.thickness; %set as constants
	md.geometry.surface=md.geometry.base+md.geometry.thickness;
    md.masstransport.spcthickness(ind)=md.geometry.thickness(ind); % thickness at groundling line stays constant throughout simulation
    

    md=transientrestart(md,2); 
	md.timestepping.start_time=2/12;

	timesteps=1/12;
	md.inversion.iscontrol=0;
	md.timestepping.final_time=50;
	md.timestepping.time_step=timesteps;
	%md.settings.output_frequency=timesteps*12;
	%md.settings.recording_frequency=timesteps*12;

	md.transient.amr_frequency=0;
	md.transient.isgroundingline=0;
	md.transient.isthermal=0;
	md.transient.ismovingfront=0;
	md.transient.requested_outputs={'default',...
		'StrainRateeffective','StrainRatexx','StrainRatexy','StrainRateyy',...
		'StrainRatexz','StrainRateyz','StrainRatezz','BasalforcingsFloatingiceMeltingRate','BasalforcingsGroundediceMeltingRate'...
        'DeviatoricStresseffective','DeviatoricStressxx','DeviatoricStressxy','DeviatoricStressyy'};

	%Parameterized melting rate
	md.basalforcings=basalforcingspico(md.basalforcings);
	md.basalforcings.basin_id=ones(md.mesh.numberofelements,1);
	md.basalforcings.num_basins=1;
	md.basalforcings.maxboxcount=5;
	md.basalforcings.isplume=1; %1 for PICOP (PICO + Plume), 0 for PICO

    %md.smb.mass_balance=zeros(md.mesh.numberofvertices,1);
    %md.basalforcings.groundedice_melting_rate=10*ones(md.mesh.numberofvertices,1).*(1-(md.mesh.y/Ly).^(1/3));
    %md.basalforcings.floatingice_melting_rate=10*ones(md.mesh.numberofvertices,1).*(1-(md.mesh.y/Ly).^(1/3));

	md.verbose=verbose('solution',true,'convergence',true);
	md.cluster=cluster;
	Smean=34.73;
	Tmean=272.53;

	md.basalforcings.farocean_salinity=[Smean Smean;timesteps md.timestepping.final_time]; %Time series
	md.basalforcings.farocean_temperature=[Tmean Tmean;timesteps md.timestepping.final_time];   %Time series
	md.miscellaneous.name='picop_transient_decoup_pert';

	md=solve(md,'tr');
	savemodel(org,md)
end
% }}}
