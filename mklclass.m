function [result, weights,variableveccell] = mklclass(xapp, yapp, xtest, ytest,kernelt,kerneloptionvect,variablevec,ranks, C)
% Example of how to use the mklsvm for  classification
%
%

%clear all
%close all
%addpath('C:\Users\yingxiali\Desktop\simplemkl');


%nbiter=1;
%ratio=0.5;
%data='ionosphere'
%C = [100];
verbose=1;

options.algo='svmclass'; % Choice of algorithm in mklsvm can be either
                         % 'svmclass' or 'svmreg'
%------------------------------------------------------
% choosing the stopping criterion
%------------------------------------------------------
options.stopvariation=0; % use variation of weights for stopping criterion 
options.stopKKT=0;       % set to 1 if you use KKTcondition for stopping criterion    
options.stopdualitygap=1; % set to 1 for using duality gap for stopping criterion

%------------------------------------------------------
% choosing the stopping criterion value
%------------------------------------------------------
options.seuildiffsigma=1e-2;        % stopping criterion for weight variation 
options.seuildiffconstraint=0.1;    % stopping criterion for KKT
options.seuildualitygap=0.01;       % stopping criterion for duality gap

%------------------------------------------------------
% Setting some numerical parameters 
%------------------------------------------------------
options.goldensearch_deltmax=1e-1; % initial precision of golden section search
options.numericalprecision=1e-8;   % numerical precision weights below this value
                                   % are set to zero 
options.lambdareg = 1e-8;          % ridge added to kernel matrix 

%------------------------------------------------------
% some algorithms paramaters
%------------------------------------------------------
options.firstbasevariable='first'; % tie breaking method for choosing the base 
                                   % variable in the reduced gradient method 
options.nbitermax=500;             % maximal number of iteration  
options.seuil=0;                   % forcing to zero weights lower than this 
options.seuilitermax=10;           % value, for iterations lower than this one 

options.miniter=0;                 % minimal number of iterations 
options.verbosesvm=0;              % verbosity of inner svm algorithm 
options.efficientkernel=1;         % use efficient storage of kernels 


%------------------------------------------------------------------------
%                   Building the kernels parameters
%------------------------------------------------------------------------
kernelt=kernelt;
kerneloptionvect=kerneloptionvect;
variablevec=variablevec;




%classcode=[1 -1];
%load([data ]);
[~,dim]=size(xapp);

%nbtrain=floor(nbdata*ratio);
%rand('state',0);


    %[xapp,yapp,xtest,ytest,indice]=CreateDataAppTest(x, y, nbtrain,classcode);
    [xapp,xtest]=normalizemeanstd(xapp,xtest);
    [kernel,kerneloptionvec,variableveccell]=CreateKernelListWithVariable(variablevec,dim,kernelt,kerneloptionvect);
    kernel=kernelt;
    kerneloptionvec=kerneloptionvect;
    variableveccell=ranks;
    [Weight,InfoKernel]=UnitTraceNormalization(xapp,kernel,kerneloptionvec,variableveccell);
    K=mklkernel(xapp,InfoKernel,Weight,options);


    
    %------------------------------------------------------------------
    % 
    %  K is a 3-D matrix, where K(:,:,i)= i-th Gram matrix 
    %
    %------------------------------------------------------------------
    % or K can be a structure with uses a more efficient way of storing
    % the gram matrices
    %
    % K = build_efficientK(K);
    
    tic
    [beta,w,b,posw,~,~] = mklsvm(K,yapp,C,options,verbose);
    toc

    Kt=mklkernel(xtest,InfoKernel,Weight,options,xapp(posw,:),beta);
    weights = Weight;
    ypred=Kt*w+b;

    %mean(sign(ypred)==ytest)
    result = ypred;

end%



