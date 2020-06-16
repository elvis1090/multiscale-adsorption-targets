function analyze_overpredictions(c0,v_total,t_total,eps,A_mem,nu,cin,rho_mem_SI,dname)
    % This function is used to analyze the amount of over-/under-predictions in
    % predicted material properties and its effect on regeneration time for a
    % given process. This analysis assumes constant v_total,t_total and cin
    % to obtain dimensional forms of the dimensionless variables
    disp(' ');
    disp('Analyzing over-/under-predictions due to scale up');

    % read data for hypothetical materials. this data file should be 
    % manually generated based on the
    % points a user would like to analyze. this data file 
    % is not generated automatically by the code. make sure to modify the
    % range accordingly if the number of hypothetical materials are changed
    % from 3
    hyp_mat_dat = readmatrix('./hypothetical_sorbent_properties.csv',...
        'Range','B2:D6');
    
    % group into meaningful variables
    K_dat = hyp_mat_dat(1,:);
    Q_bat_dat = hyp_mat_dat(2,:);
    Q_semicont_dat = hyp_mat_dat(3,:);
    r_dat = hyp_mat_dat(4,:);
    nbv_dat = hyp_mat_dat(5,:);
    
    % use cin to obtain values with dimensions for use with our functions
    K_wdim = K_dat / cin;
    Q_bat_wdim = Q_bat_dat * cin;
    Q_semicont_wdim = Q_semicont_dat * cin;
    cout_wdim = cin./r_dat;
    v_mem_wdim = v_total./nbv_dat;
    
    % Calculate module masses and length
    m_mem = v_mem_wdim*rho_mem_SI;
    mod_len = v_mem_wdim/A_mem;
    
    % Get number of hypothetical materials
    n_hyp_mats = length(K_dat);
    
    % define an anonymous function to pass constant data to the
    % semicont_breakthrough function
    breakthrough = @(K,Q,clim,module_len,module_mass) semicont_breakthrough(K,Q,10,'suff',cin,v_total,t_total,eps,A_mem,nu,module_len,c0,clim,module_mass,rho_mem_SI,dname);
                                                      
    % Iterate over hypothetical materials and calculate breakthrough times
    % and bed volumes
    for i = 1:n_hyp_mats
        disp(['Analyzing hypothetical material ',num2str(i)])
                
        % breakthrough times
        disp(['Breakthrough time calculations'])
        disp('Breakthrough time from batch upscaling')
        bt_batch_pred = breakthrough(K_wdim(i),Q_bat_wdim(i),cout_wdim(i),mod_len(i),m_mem(i));
        disp('Breakthrough time from semi-continuous upscaling')
        bt_semicont_pred = breakthrough(K_wdim(i),Q_semicont_wdim(i),cout_wdim(i),mod_len(i),m_mem(i));
        
        
        % calculate and display underprediction
        underpred = (bt_semicont_pred-bt_batch_pred)/bt_semicont_pred * 100;
        disp(['Underprediction = ',num2str(underpred),'%'])
        disp('-------')
        
        % bed volumes
        disp(['Bed volume calculations'])
        bat_nbv_dat = batch_nbv(K_dat(i), Q_bat_dat(i), 0.3, r_dat(i));
        disp(['Kbar=',num2str(K_dat(i)),' Qbar=',num2str(Q_bat_dat(i)),' r=',num2str(r_dat(i))]);
        disp(['N_BV Batch = ',num2str(bat_nbv_dat)]);
        cont_nbv_dat = cont_nbv(K_dat(i), Q_bat_dat(i), 0.3);
        disp(['Kbar=',num2str(K_dat(i)),' Qbar=',num2str(Q_bat_dat(i))]);
        disp(['N_BV Semi-continuous=',num2str(cont_nbv_dat)]);
        nbv_overpredict = (bat_nbv_dat-cont_nbv_dat)/bat_nbv_dat * 100;
        display(['Overprediction (w.r.t) batch = ',num2str(nbv_overpredict),'%']);
        disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    end
end