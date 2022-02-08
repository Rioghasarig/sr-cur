classdef lusol_obj < handle
  %LUSOL_OBJ  access to LUSOL for computations on a sparse LU factorization.
  %
  % Optional parameters:
  %  options = lusol_obj.luset();
  %  options = lusol_obj.luset('pivot','TRP');
  %
  % Initialization:
  %  mylu = lusol_obj(A);
  %
  % Initialization specifying options:
  %  mylu = lusol_obj(A,options);
  %
  % Note that initialization performs the factorization.  Subsequent
  % factorizations with the same object may be performed with the factorize
  % method.
  %
  % Factorize:
  %  [inform nsing depcol] = mylu.factorize(A);
  %  [inform nsing depcol] = mylu.factorize(A,options);
  %
  % There are several operations that may be performed with a lusol object.
  %
  % Usage:
  %  y = mylu.mulA(x);   % y = A*x
  %  y = mylu.mulAt(x);  % y = A'*x
  %  x = mylu.solveA(b); % x = A\b
  %  x = mylu.solveAt(b); % x = A'\b
  %
  % Update:
  %  inform = mylu.repcol(v,j); % replace column j of A with vector v
  %
  % See also:
  %  lusol
  %  lusol_obj.luset
  %  lusol_obj.factorize
  %  lusol_obj.stats
  %  lusol_obj.L0
  %  lusol_obj.U
  %  lusol_obj.p
  %  lusol_obj.q
  %  lusol_obj.mul
  %  lusol_obj.mulA
  %  lusol_obj.mulAt
  %  lusol_obj.mulL
  %  lusol_obj.mulLt
  %  lusol_obj.mulU
  %  lusol_obj.mulUt
  %  lusol_obj.solve
  %  lusol_obj.solveA
  %  lusol_obj.solveAt
  %  lusol_obj.solveL
  %  lusol_obj.solveLt
  %  lusol_obj.solveU
  %  lusol_obj.solveUt
  %  lusol_obj.repcol
  %  lusol_obj.reprow
  %  lusol_obj.addcol
  %  lusol_obj.addrow
  %  lusol_obj.delcol
  %  lusol_obj.delrow
  %  lusol_obj.r1mod
  %

  properties (Access=public)

    % object parameters

    nzinit = 0; % initial number of non zeros
    nzinit2 = 0; 
    
    % lusol input parameters
    rank = 0;
    maxcol = 0; % max num cols searched for piv element (5)
    pivot = 0; % pivoting method [(0=TPP),1=TRP,2=TCP,3=TSP]
    keepLU = 0; % keep the nonzeros, if 0, permutations are computed (1)
    Ltol1 = 0; % max Lij allowed, default depends on luparm(6)
    Ltol2 = 0; % max Lij allowed during updates
    small = 0; % absolute tolerance for treating reals as zero (eps^0.8)
    Utol1 = 0; % absolute tol for flagging small diags of U (eps^0.67)
    Utol2 = 0; % rel tol for flagging small diags of U (eps^0.67)
    Uspace = 0; % (3.0)
    dens1 = 0; % (0.3)
    dens2 = 0; % (0.5)

    % lusol parameter vectors

    luparm_ptr = 0; % vector of integer parameters (input and output)
    parmlu_ptr = 0; % vector of double parameters (input and output)

    % scalars

    m_ptr = 0; % number of rows
    n_ptr = 0; % number of columns

    nelem_ptr = 0; % number of elements in original matrix
    nzmax_ptr = 0; % maximum storage allocated
    rank_ptr = 0; % rank of the factorization 
    % vectors of lenth nzmax

    a_ptr = 0; % main storage array
    indc_ptr = 0; % row indices
    indr_ptr = 0; % column indices

    % vectors of length m

    p_ptr = 0; % row permutation
    ap_ptr = 0;
    lenr_ptr = 0;
    locr_ptr = 0;
    

    iqloc_ptr = 0;
    ipinv_ptr = 0;
    ap = []
    % vectors of length n

    q_ptr = 0; % column permutation
    aq_ptr = 0; 
    lenc_ptr = 0;
    locc_ptr = 0;
    
    iploc_ptr = 0;
    iqinv_ptr = 0;
    aq = [];
    
    depcol_lx = 0; % logical index indicating dependent columns
    int_class = 'int64'; % integer class used for integer arrays
    int_ptr_class = 'int64Ptr'; % integer class for libpointers
    
    A = [];
    L21 = [];
    U12t = [];
    L11 = [];
    U11 = []; 
    U12 = [];
    
    lenlv_ptr = 0;
    li_ptr = 0;
    lj_ptr = 0;
    lv_ptr = 0;
    
    m = 0;
    n = 0;
    
    O = [];
    S = [];
    A11 = [];
    A22 = [];
    A12 = [];
    A21 = [];
    B21 = []; 
    
    lops = [];
    li = [];
    lj = [];
    lv = []; 
    
    Qc = [];
    Bk = [];
    B = []; 
    Qr = []; 
    
    Rc = [];
    Rr = []; 
  end

  methods (Static)

    function options = luset(varargin)
      %LUSET  process input to create lusol options structure
      %
      % This method uses Matlab's inputParser class to handle lusol options.
      % It handles empty input, input structures, and key-value lists just like
      % many other Matlab programs.
      %
      % To obtain a structure with default settings use:
      %
      %  options = lusol.luset();
      %
      % To create an options structre that will use threshold rook pivoting,
      % use:
      %
      %  options = lusol.luset('pivot','TRP')
      %
      % The best reference for the lusol parameters is currently the comments
      % for the lu1fac subroutine in lusol1.f.
      %
      % |--------+----------+----------------------------------------------------|
      % | param  |  default | description                                        |
      % |--------+----------+----------------------------------------------------|
      %
      % lusol_obj options
      % |--------+----------+----------------------------------------------------|
      % | nzinit |        0 | minimum length for storage arrays                  |
      % |--------+----------+----------------------------------------------------|
      %
      % LUSOL integer parameters
      % |--------+----------+----------------------------------------------------|
      % | maxcol |        5 | max num cols searched for piv element              |
      % | pivot  |    'TPP' | pivoting method {'TPP','TRP','TCP','TSP'}          |
      % | keepLU |        1 | keep the nonzeros, if 0, permutations are computed |
      % |--------+----------+----------------------------------------------------|
      %
      % LUSOL real parameters
      % |--------+----------+----------------------------------------------------|
      % | Ltol1  |     10.0 | max Lij allowed, default depends on pivot method   |
      % | Ltol2  |     10.0 | max Lij allowed during updates                     |
      % | small  |  eps^0.8 | absolute tolerance for treating reals as zero      |
      % | Utol1  | eps^0.67 | absolute tol for flagging small diags of U         |
      % | Utol2  | eps^0.67 | rel tol for flagging small diags of U              |
      % | Uspace |      3.0 |                                                    |
      % | dens1  |      0.3 |                                                    |
      % | dens2  |      0.5 |                                                    |
      % |--------+----------+----------------------------------------------------|

      % get the input parser
      in_parse = inputParser;

      % storage parameters
      in_parse.addParamValue('nzinit',0,@(x) x>=0);
      in_parse.addParamValue('nzinit2',0,@(x) x>=0); 
      
      % lusol integer parameters
      in_parse.addParamValue('maxcol',5,@(x) x>=0);
      in_parse.addParamValue('pivot','TPP',@(x) ismember(x,{'TPP','TRP','TCP','TSP'}));
      in_parse.addParamValue('keepLU',1,@(x) ismember(x,[0 1]));

      % lusol real parameters
      in_parse.addParamValue('Ltol1',10.0,@(x) x>=0.0);
      in_parse.addParamValue('Ltol2',10.0,@(x) x>=0.0);
      in_parse.addParamValue('small',eps^0.8,@(x) x>=0.0);
      in_parse.addParamValue('Utol1',eps^0.67,@(x) x>=0.0);
      in_parse.addParamValue('Utol2',eps,@(x) x>=0.0);
      in_parse.addParamValue('Uspace',3.0,@(x) x>=0.0);
      in_parse.addParamValue('dens1',0.3,@(x) x>=0.0);
      in_parse.addParamValue('dens2',0.5,@(x) x>=0.0);
       

      in_parse.addParamValue('rank',1, @(x) x>=0); 
      % parse the input
      in_parse.parse(varargin{:});

      % obtain the output
      options = in_parse.Results;

    end

    function r = get_range(p,i,j)
      %GET_RANGE  return a subset array from a libpointer array

      % get the datatype of the libpointer object
      pdatatype = p.DataType;
      % create a new pointer to p[i-1]
      pp = p + (i-1);
      % change the size of pp
      pp.setdatatype(pdatatype,j-i+1);
      % get the return array
      r = pp.Value;

    end

    function v = get_value(p,i)
      %GET_VALUE  return a single value from a libpointer array

      v = lusol_obj.get_range(p,i,i);

    end

    function load_library
      %LOAD_LIBRARY  Load the clusol shared library

      % do nothing if the library is already loaded
      if libisloaded('libclusol')
        return;
      end

      switch lower(computer)
        case 'glnxa64'
          % load clusol with linux prototype file
          loadlibrary('libclusol',@libclusol_proto_glnxa64);
        case 'maci64'
          % load clusol with mac prototype file
          loadlibrary('libclusol',@libclusol_proto_maci64);
        otherwise
          % the interface has not been implemented for other systems
          error('lusol_obj:load_library', ...
            'clusol library not implemented for this architecture')
      end

    end

    function unload_library
      %UNLOAD_LIBRARY  Unload the clusol shared library

      if libisloaded('libclusol')
        unloadlibrary('libclusol');
      end

    end

    function y = vector_process(x,xn)
      %VECTOR_PROCESS check and prepare vector for use in update routines
      %
      % Usage:
      %  y = lusol_obj.vector_process(x,xn);
      %
      % Input:
      %  x = input vector
      %  xn = expected length of x
      %
      % Output:
      %  y = output vector (double, full, and columnar)
      %
      % Error:
      %  The function will throw an error if x is not a vector of length xn
      %

      % check input vector
      if ~isvector(x) || length(x) ~= xn
        error('lusol_obj:vector_process','input is not a vector of correct length.');
      end

      % orient and densify
      y = double(full(x(:)));

    end

  end

  methods (Access=private)

    function parse_options(obj,varargin)
      %PARSE_OPTIONS  process the options structure to set parameters.

      % use luset to parse user input
      options = obj.luset(varargin{:});

      % storage parameters
      obj.nzinit = options.nzinit;
      obj.nzinit2 = options.nzinit2;
      
      % lusol integer parameters
      obj.maxcol = options.maxcol;
      obj.keepLU = options.keepLU;
      if(obj.rank == 0)
        obj.rank = options.rank;
      end
      % lusol double parameters
      obj.Ltol1 = options.Ltol1;
      obj.Ltol2 = options.Ltol2;
      obj.small = options.small;
      obj.Utol1 = options.Utol1;
      obj.Utol2 = options.Utol2;
      obj.Uspace = options.Uspace;
      obj.dens1 = options.dens1;
      obj.dens2 = options.dens2;

      % set the pivoting strategy
      switch options.pivot
        case 'TPP'
          obj.pivot = 0;
        case 'TRP'
          obj.pivot = 1;
        case 'TCP'
          obj.pivot = 2;
        case 'TSP'
          obj.pivot = 3;
        otherwise
          error('lusol:set_options','Unkown pivot strategy.')
      end
    end

    function set_options(obj)
      %SET_OPTIONS  allocate and assign parameters to LUSOL arrays
      %
      % LUSOL stores input and output scalar parameters in two vectors:
      %
      %  luparm is an integer array of length 30
      %  parmlu is a double array of length 30
      %
      % this method sets the LUSOL input parameters in the correct location
      % for the fortran calls.
      %

      % allocate parameter vectors
      luparm = zeros(30,1,obj.int_class);
      parmlu = zeros(30,1,'double');

      % set parameter values
      luparm(2) = cast(-1,obj.int_class);
      luparm(3) = cast(obj.maxcol,obj.int_class);
      luparm(6) = cast(obj.pivot,obj.int_class);
      luparm(8) = cast(obj.keepLU,obj.int_class);
      parmlu(1) = double(obj.Ltol1);
      parmlu(2) = double(obj.Ltol2);
      parmlu(3) = double(obj.small);
      parmlu(4) = double(obj.Utol1);
      parmlu(5) = double(obj.Utol2);
      parmlu(6) = double(obj.Uspace);
      parmlu(7) = double(obj.dens1);
      parmlu(8) = double(obj.dens2);

      % allocate and initialize pointers
      obj.luparm_ptr = libpointer(obj.int_ptr_class,luparm);
      obj.parmlu_ptr = libpointer('doublePtr',parmlu);

    end

    function allocate_and_copy(obj,A)
      %ALLOCATE_AND_COPY  allocate LUSOL storage arrays
      %
      % LUSOL operates on many arrays.  This method allocates all of them
      % to an appropriate size.
      %

      % get information about A
      m = size(A,1);
      n = size(A,2);

      nelem = nnz(A);

      % set storage sizes
     
      nzmax = max([2*nelem 10*m 10*n 10000 obj.nzinit]);
      rank_ = max(obj.rank, 1); 
      % vectors of length nzmax
      a = zeros(nzmax,1);
      indc = zeros(nzmax,1,obj.int_class);
      indr = zeros(nzmax,1,obj.int_class);
      
      % extract data from A for use in LUSOL
      [indc_tmp indr_tmp a_tmp] = find(A);
      indc(1:nelem) = cast(indc_tmp,obj.int_class);
      indr(1:nelem) = cast(indr_tmp,obj.int_class);
      a(1:nelem) = a_tmp;
      % vectors of length m
      p = zeros(m,1);
      ap = zeros(m,1); 

      lenr = zeros(m,1);
      locr = zeros(m,1);
      

      iqloc = zeros(m,1);
      ipinv = zeros(m,1);
      % vectors of length n
      q = zeros(n,1);
      aq  = zeros(n,1); 

      lenc = zeros(n,1);
      locc = zeros(n,1);
      
      iploc = zeros(n,1);
      iqinv = zeros(n,1);

      li = zeros(3*rank_,1);
      lj = zeros(3*rank_,1);
      lv = zeros(3*rank_,1);
      
      %-- allocate and initialize libpointer "arrays" --%
      % integer scalars
      obj.m_ptr = libpointer(obj.int_ptr_class,m);
      obj.n_ptr = libpointer(obj.int_ptr_class,n);
      obj.nelem_ptr = libpointer(obj.int_ptr_class,nelem);
      obj.nzmax_ptr = libpointer(obj.int_ptr_class,nzmax);
      obj.rank_ptr = libpointer(obj.int_ptr_class, rank_); 
      % vectors of length nzmax
      obj.a_ptr = libpointer('doublePtr',a);
      obj.indc_ptr = libpointer(obj.int_ptr_class,indc);
      obj.indr_ptr = libpointer(obj.int_ptr_class,indr);


      % vectors of length m
      obj.p_ptr = libpointer(obj.int_ptr_class,p);
      obj.ap_ptr = libpointer(obj.int_ptr_class,ap); 
      
      obj.lenr_ptr = libpointer(obj.int_ptr_class,lenr);
      obj.locr_ptr = libpointer(obj.int_ptr_class,locr);
      

      obj.iqloc_ptr = libpointer(obj.int_ptr_class,iqloc);
      obj.ipinv_ptr = libpointer(obj.int_ptr_class,ipinv);
      % vectors of length n
      obj.q_ptr = libpointer(obj.int_ptr_class,q);
      obj.aq_ptr = libpointer(obj.int_ptr_class,aq); 

      obj.lenc_ptr = libpointer(obj.int_ptr_class,lenc);
      obj.locc_ptr = libpointer(obj.int_ptr_class,locc);
      
      obj.iploc_ptr = libpointer(obj.int_ptr_class,iploc);
      obj.iqinv_ptr = libpointer(obj.int_ptr_class,iqinv);


      obj.lenlv_ptr = libpointer(obj.int_ptr_class, [0,0,0]);
      obj.li_ptr = libpointer(obj.int_ptr_class, li);
      obj.lj_ptr = libpointer(obj.int_ptr_class, lj);
      obj.lv_ptr = libpointer('doublePtr', lv);
      
    end

    function update_check(obj)
      %UPDATE_CHECK  throw an error if this method is called after updates
      %
      % Some methods should not be called after updates to a factorization.
      % This method checks if any updated have occuered and throws and
      % error if this is the case.
      %

      nupdat = lusol_obj.get_value(obj.luparm_ptr,15);
      if nupdat > 0
        error('lusol:post_update_call_error', ...
              'nsing and depcol cannot be called after updates.')
      end

    end

    function [x inform resid] = clu6sol(obj,b,mode)
      %CLU6SOL  call lu6sol to perform various solves with L and U factors.
      %
      % Performs data processing and direct call to clusol function for
      % solves.  Right hand side must be a vector.
      %
      % This is a private method and should only be used by class methods.
      % Users should call the various public interface methods.
      %
      % Usage:
      %  x = obj.clu6sol(b,mode)
      %
      % Input:
      %  b = right hand side vector
      %  mode = solution mode (see table below)
      %
      % Output:
      %  x = solution vector
      %  inform = status flag
      %  resid = 1-norm of residual
      %
      % Modes:
      %  1    x  solves   L x = b
      %  2    x  solves   L'x = b
      %  3    x  solves   U x = b
      %  4    x  solves   U'x = b
      %  5    x  solves   A x = b (default)
      %  6    x  solves   A'x = b
      %
      % inform flags:
      %  0 = successful solve
      %  1 = if U is singular, and residual is non-zero
      %

      % handle function options
      if nargin < 3
        mode = 5;
      end
      % make sure b is a vector
      if ~isvector(b)
        error('lusol:solve','b must be a vector.')
      end
      % orient b vector
      b = double(full(b(:)));
      lenb = length(b);
      % get size
      [m n] = obj.size();
      % allocate v and w vectors
      v = zeros(m,1,'double');
      w = zeros(n,1,'double');
      switch mode
        case 1
          if lenb ~= m, error('lusol:solve','b has incorrect size.'); end
          v = b;
        case 2
          if lenb ~= m, error('lusol:solve','b has incorrect size.'); end
          v = b;
        case 3
          if lenb ~= m, error('lusol:solve','b has incorrect size.'); end
          v = b;
        case 4
          if lenb ~= n, error('lusol:solve','b has incorrect size.'); end
          w = b;
        case 5
          if lenb ~= m, error('lusol:solve','b has incorrect size.'); end
          v = b;
        case 6
          if lenb ~= n, error('lusol:solve','b has incorrect size.'); end
          w = b;
        otherwise
          error('lusol:solve','unrecognized mode.')
      end

      % set up local libpointers for function call
      v_ptr = libpointer('doublePtr',v);
      w_ptr = libpointer('doublePtr',w);
      mode_ptr = libpointer(obj.int_ptr_class,mode);
      ret_inform_ptr = libpointer(obj.int_ptr_class,0);

      % call lusol routine
      calllib('libclusol','clu6sol', ...
        mode_ptr, ...
        obj.m_ptr, ...
        obj.n_ptr, ...
        v_ptr, ...
        w_ptr, ...
        obj.nzmax_ptr, ...
        obj.luparm_ptr, ...
        obj.parmlu_ptr, ...
        obj.a_ptr, ...
        obj.indc_ptr, ...
        obj.indr_ptr, ...
        obj.p_ptr, ...
        obj.q_ptr, ...
        obj.lenc_ptr, ...
        obj.lenr_ptr, ...
        obj.locc_ptr, ...
        obj.locr_ptr, ...
        ret_inform_ptr);

      switch mode
        case 1
          x = v_ptr.Value;
        case 2
          x = v_ptr.Value;
        case 3
          x = w_ptr.Value;
        case 4
          x = v_ptr.Value;
        case 5
          x = w_ptr.Value;
        case 6
          x = v_ptr.Value;
      end

      inform = double(lusol_obj.get_value(obj.luparm_ptr,10));
      resid = double(lusol_obj.get_value(obj.parmlu_ptr,20));
    end

    function y = clu6mul(obj,x,mode)
      %CLU6MUL  call LUSOL to perform various multiplies with L and U factors.
      %
      % Performs data processing and direct call to clusol function to compute
      % a matrix vector multiply with the desired factor.
      %
      % This is a private method and should only be used by class methods.
      % Users should call the various public interface methods.
      %
      % Usage:
      %  y = obj.clu6mul(x,mode)
      %
      % Input:
      %  x = vector to multiply
      %  mode = multiply mode, see table below
      %
      % Output:
      %  y = product of desired factor and x
      %
      % Modes:
      %  1    y = L x
      %  2    y = L'x
      %  3    y = U x
      %  4    y = U'x
      %  5    y = A x (default)
      %  6    y = A'x
      %

      % handle optional input
      if nargin < 3
        mode = 5;
      end

      % check if x is a vector
      if ~isvector(x)
        error('lusol:mul','x must be a vector.');
      end

      % orient x vector
      x = double(full(x(:)));
      lenx = length(x);

      % get matrix size
      [m n] = obj.size();

      % set up temporary vectors
      v = zeros(m,1);
      w = zeros(n,1);

      switch mode
        case 1
          if lenx ~= m, error('lusol:mul','x has incorrect size.'); end
          v = x;
        case 2
          if lenx ~= m, error('lusol:mul','x has incorrect size.'); end
          v = x;
        case 3
          if lenx ~= n, error('lusol:mul','x has incorrect size.'); end
          w = x;
        case 4
          if lenx ~= m, error('lusol:mul','x has incorrect size.'); end
          v = x;
        case 5
          if lenx ~= n, error('lusol:mul','x has incorrect size.'); end
          w = x;
        case 6
          if lenx ~= m, error('lusol:mul','x has incorrect size.'); end
          v = x;
        otherwise
          error('lusol:mul','unrecognized mode.')
      end

      % set up local libpointers
      mode_ptr = libpointer(obj.int_ptr_class,mode);
      v_ptr = libpointer('doublePtr',v);
      w_ptr = libpointer('doublePtr',w);

      % call clu6mul
      calllib('libclusol','clu6mul', ...
        mode_ptr, ...
        obj.m_ptr, ...
        obj.n_ptr, ...
        v_ptr, ...
        w_ptr, ...
        obj.nzmax_ptr, ...
        obj.luparm_ptr, ...
        obj.parmlu_ptr, ...
        obj.a_ptr, ...
        obj.indc_ptr, ...
        obj.indr_ptr, ...
        obj.p_ptr, ...
        obj.q_ptr, ...
        obj.lenc_ptr, ...
        obj.lenr_ptr, ...
        obj.locc_ptr, ...
        obj.locr_ptr);

      switch mode
        case 1
          y = v_ptr.Value;
        case 2
          y = v_ptr.Value;
        case 3
          y = v_ptr.Value;
        case 4
          y = w_ptr.Value;
        case 5
          y = v_ptr.Value;
        case 6
          y = w_ptr.Value;
      end

    end

  end

  methods

    % constructor and main factorize method

    function obj = lusol_obj(A,varargin)
      %LUSOL_OBJ  constructor for lusol object, factorize A
      %
      % Creates lusol object and factorizes A.
      %
      % Example:
      %   mylu = lusol(A);
      %
      % See also: LUSOL_OBJ.FACTORIZE
      %

      % load the shared library
      obj.load_library;

      % factorize the matrix
      obj.factorize(A,varargin{:});
      [m,n] = size(A); 
      obj.m = m;
      obj.n = n;
      obj.A = A;
      [L,U] = obj.set_data();
      nrank = obj.rank;

      % ===================
      obj.ap = obj.ap_ptr.Value;
      obj.aq = obj.aq_ptr.Value;
      
      obj.A11 = obj.A(obj.ap(1:nrank),obj.aq(1:nrank));
      obj.A12 = obj.A(obj.ap(1:nrank),obj.aq(nrank+1:end));
      obj.A21 = obj.A(obj.ap(nrank+1:end),obj.aq(1:nrank));
      obj.A22 = obj.A(obj.ap(nrank+1:end),obj.aq(nrank+1:end));
      
      obj.L21 = L(nrank+1:end,1:nrank);       
      obj.U12t = U(1:nrank,nrank+1:end)';
      
      obj.O = randn(20,obj.m-nrank);
      OLU = obj.O*obj.L21;
      OLU = OLU*obj.U12t';
      obj.S = obj.O*obj.A(obj.ap(nrank+1:end),obj.aq(nrank+1:end)) - ...
                OLU;
            
            
    end

    function r = factorize(obj,A,varargin)
      %FACTORIZE  perform lu factorization on A
      %
      % This method tells LUSOL to perform an LU factorization on A.
      %
      % Usage (after mylu object is initialized):
      %  [inform nsing depcol] = mylu.factorize(A)
      %  [inform nsing depcol] = mylu.factorize(A,options)
      %
      % Input:
      %  A = matrix to factorize
      %  options = options structure (optional)
      %
      % Output:
      %  inform = status flag
      %  nsing = esimate of the number of singularities
      %  depcol = logical index of dependent columns
      %
      % inform code:
      %  0 if the LU factors were obtained successfully.
      %  1 if U appears to be singular, as judged by lu6chk.
      %  3 if some index pair indc(l), indr(l) lies outside
      %    the matrix dimensions 1:m , 1:n.
      %  4 if some index pair indc(l), indr(l) duplicates
      %    another such pair.
      %  7 if the arrays a, indc, indr were not large enough.
      %    Their length "lena" should be increased to at least
      %    the value "minlen" given in luparm(13).
      %  8 if there was some other fatal error.  (Shouldn't happen!)
      %  9 if no diagonal pivot could be found with TSP or TDP.
      %    The matrix must not be sufficiently definite
      %    or quasi-definite.
      %

      % parse and set options
      obj.parse_options(varargin{:});
      obj.set_options();

      % allocate arrays and copy data from matrix A
      obj.allocate_and_copy(A);

      % get matrix size
      [m n] = obj.size();
      
      % temporary storage
      %cols_ptr = libpointer(obj.int_ptr_class,zeros(n,1));
      %markc_ptr = libpointer(obj.int_ptr_class,zeros(n,1));
      %markr_ptr = libpointer(obj.int_ptr_class,zeros(m,1));
      w_ptr = libpointer('doublePtr',zeros(n,1));

      % run lusol
      ret_inform_ptr = libpointer(obj.int_ptr_class,0);
      calllib('libclusol','clu1fac', ...
        obj.m_ptr, ...
        obj.n_ptr, ...
        obj.nelem_ptr, ...
        obj.nzmax_ptr, ...
        obj.ap_ptr, ...
        obj.aq_ptr, ...
        obj.rank_ptr, ...
        obj.luparm_ptr, ...
        obj.parmlu_ptr, ...
        obj.a_ptr, ...
        obj.indc_ptr, ...
        obj.indr_ptr, ...
        obj.p_ptr, ...
        obj.q_ptr, ...
        obj.lenc_ptr, ...
        obj.lenr_ptr, ...
        obj.locc_ptr, ...
        obj.locr_ptr, ...
        obj.iploc_ptr, ...
        obj.iqloc_ptr, ...
        obj.ipinv_ptr, ...
        obj.iqinv_ptr, ...
        w_ptr, ...
        ret_inform_ptr);

      % error checking
      ret_inform = ret_inform_ptr.Value;
      switch ret_inform
        case 0
          % ok, LU factors obtained
        case 1
          % ok, LU factors obtained, rank deficient
        case 3
          % error, some index pair indc(l), indr(l) lies outside
          % the matrix dimensions 1:m , 1:n.  Should not happen, because
          % matlab controls the input to lu1fac.
          err = MException('lusol:factorize', ...
                           'LUSOL reports improper input. inform = %d',ret_inform);
          throw(err);
        case 4
          % error, some index pair indc(l), indr(l) duplicates
          % another such pair.  Should not happen, because
          % matlab controls the input to lu1fac.
          err = MException('lusol:factorize', ...
                           'LUSOL reports improper input. inform = %d',ret_inform);
          throw(err);
        case 7
          % error, not enough storage.  User needs to increase nzinit
          % parameter
          err = MException('lusol:factorize', ...
                           ['LUSOL needs more storage. ' ...
                            'Increase the nzinit parameter.  inform = %d'], ...
                           ret_inform);
          throw(err);
        case 8
          % some other fatal error
          err = MException('lusol:factorize', ...
                           'LUSOL other fatal error. inform = %d',ret_inform);
          throw(err);
        case 9
          % error, no diagonal pivot could be found with TSP or TDP.
          % The matrix must not be sufficiently definite or quasi-definite
          err = MException('lusol:factorize', ...
                           ['LUSOL no diagonal pivot could be found with ' ...
                            'TSP or TDP. inform = %d'], ...
                           ret_inform);
          throw(err);
      end

      % compute logical index of dependent columns
      obj.depcol_lx = (w_ptr.Value <= 0.0);

      % user requests inform flag
      if nargout > 0
        inform = obj.inform();
      end

      % user requests number of singularities
      if nargout > 1
        nsing = obj.nsing();
      end

      % user requests dependent column indicator
      if nargout > 2
        depcol = obj.depcol_lx;
      end

    end
    
    function [L,U] = set_data(obj)
      nrank = obj.rank; 
      
      % Extract matrices
      L = obj.L0();
      U = obj.getU11();

      
      U11 = U(1:nrank,1:nrank); 
      obj.U11 = U11;
      [uj,~,uv] = find(U11'); 
      lenu = length(uv); 
      
      L11 = L(1:nrank,1:nrank);
      obj.L11 = L11;
      L11 = L11 - speye(nrank);
      [i,j,v] = find(L11); 
      lenl = length(v);
      
          
      nelem = lenu + lenl;
      nzmax = max([2*nelem 10*obj.m 10*obj.n 5*10^6 obj.nzinit2]);
      a = zeros(nzmax,1);
      indc = zeros(nzmax,1);
      indr = zeros(nzmax,1); 
      lena = nzmax;   
      luparm = obj.luparm_ptr.Value;
      
      % Resize factorization
      m = nrank;
      n = nrank; 
      p = 1:nrank;
      q = 1:nrank;
      
      %Extract L and U

      if(lena <= lenu + lenl)
          error('lusol:lusol_obj', "Insufficient storage");
      end
      
      % Set U data
      a(1:lenu) = uv;
      indr(1:lenu) = uj;          
      lenr = full(sum(U11 ~= 0,2))';
      locr = cumsum([1,lenr(1:nrank-1)]);
      
            
      luparm(22) = lenu; % set lenU0 to lenu
      luparm(24) = lenu; % set lenU to lenu
      luparm(25) = lenu; % set lrow to lenu
      
      %Set L data
      a(lena-lenl+1:end) = flip(-v);
      indc(lena-lenl+1:end) = flip(i);
      indr(lena-lenl+1:end) = flip(j);
      
      luparm(20) = 0; % set numL0 to 0
      luparm(21) = 0; % set lenL0 to 0
      luparm(23) = lenl; %set lenL to lenl
      
      % Set pointers
      obj.m_ptr = libpointer(obj.int_ptr_class, m); 
      obj.n_ptr = libpointer(obj.int_ptr_class, n);

      obj.p_ptr = libpointer(obj.int_ptr_class,p); 
      obj.q_ptr = libpointer(obj.int_ptr_class,q);
      
      
      obj.a_ptr = libpointer('doublePtr',a);
      obj.indc_ptr = libpointer(obj.int_ptr_class,indc);
      obj.indr_ptr = libpointer(obj.int_ptr_class,indr);
      
      obj.locr_ptr = libpointer(obj.int_ptr_class, locr);
      obj.lenr_ptr = libpointer(obj.int_ptr_class, lenr);      

      obj.luparm_ptr = libpointer(obj.int_ptr_class,luparm); 
      obj.nzmax_ptr = libpointer(obj.int_ptr_class,nzmax); 
    end
    
    function [inform nsing depcol] = refactor(obj,varargin)
      
      nrank = obj.rank;
      obj.factorize(obj.A11, varargin{:});
      [L11,U11] = obj.set_data();
      if( obj.stats.nrank ~= nrank)
          error('lusol:refactor','A11 not full rank');
      end
        
      np = obj.ap_ptr.Value; 
      nq = obj.aq_ptr.Value;
        
      obj.ap(1:nrank) = obj.ap(np); 
      obj.aq(1:nrank) = obj.aq(nq); 
        
      obj.A11 = obj.A11(np,nq);
      obj.A12 = obj.A12(np,:);
      obj.A21 = obj.A21(:,nq); 
        
      obj.L11 = L11;
      obj.U11 = U11; 
        
      obj.U12t = (obj.L11 \ obj.A12)';
      obj.L21 = obj.A21 / obj.U11; 
        
      %obj.O = randn(20,obj.m-nrank);
      OLU = obj.O*obj.L21;
      OLU = OLU*obj.U12t';
      obj.S = obj.O*obj.A(obj.ap(nrank+1:end),obj.aq(nrank+1:end)) - ...
                OLU;
                  
    end

    % methods to collect information about matrix and factorization

    function [m n] = size(obj)
      %SIZE  get size of factorized matrix
      m = double(obj.m_ptr.Value);
      n = double(obj.n_ptr.Value);
    end

    function s = stats(obj)
      %STATS  return LUSOL stats structure
      %
      % This method builds a Matlab struct containing LUSOL output
      % parameters.
      %
      % LUSOL output parameters:
      %
      % inform   Return code from last call to any LU routine.
      % nsing    No. of singularities marked in the
      %          output array w(*).
      % jsing    Column index of last singularity.
      % minlen   Minimum recommended value for  lena.
      % maxlen   ?
      % nupdat   No. of updates performed by the lu8 routines.
      % nrank    No. of nonempty rows of U.
      % ndens1   No. of columns remaining when the density of
      %          the matrix being factorized reached dens1.
      % ndens2   No. of columns remaining when the density of
      %          the matrix being factorized reached dens2.
      % jumin    The column index associated with DUmin.
      % numL0    No. of columns in initial  L.
      % lenL0    Size of initial  L  (no. of nonzeros).
      % lenU0    Size of initial  U.
      % lenL     Size of current  L.
      % lenU     Size of current  U.
      % lrow     Length of row file.
      % ncp      No. of compressions of LU data structures.
      % mersum   lu1fac: sum of Markowitz merit counts.
      % nUtri    lu1fac: triangular rows in U.
      % nLtri    lu1fac: triangular rows in L.
      % Amax     Maximum element in  A.
      % Lmax     Maximum multiplier in cu

      % Umax     Maximum element in current  U.
      % DUmax    Maximum diagonal in  U.
      % DUmin    Minimum diagonal in  U.
      % Akmax    Maximum element generated at any stage
      %          during TCP factorization.
      % growth   TPP: Umax/Amax    TRP, TCP, TSP: Akmax/Amax.
      % resid    lu6sol: residual after solve with U or U'.
      %

      luparm = double(obj.luparm_ptr.Value);
      parmlu = double(obj.parmlu_ptr.Value);

      s.inform = luparm(10);
      s.nsing = luparm(11);
      s.jsing = luparm(12);
      s.minlen = luparm(13);
      s.maxlen = luparm(14);
      s.nupdat = luparm(15);
      s.nrank = luparm(16);
      s.nrank1 = luparm(17);
      s.ndens2 = luparm(18);
      s.jumin = luparm(19);
      s.numL0 = luparm(20);
      s.lenL0 = luparm(21);
      s.lenU0 = luparm(22);
      s.lenL = luparm(23);
      s.lenU = luparm(24);
      s.lrow = luparm(25);
      s.ncp = luparm(26);
      s.mersum = luparm(27);
      s.nUtri = luparm(28);
      s.nLtri = luparm(29);
      s.Amax = parmlu(10);
      s.Lmax = parmlu(11);
      s.Umax = parmlu(12);
      s.DUmax = parmlu(13);
      s.DUmin = parmlu(14);
      s.Akmax = parmlu(15);
      s.growth = parmlu(16);
      s.resid = parmlu(20);

    end

    function info = inform(obj)
      %INFORM  return code from last call to LUSOL routines

      %info = double(obj.luparm_ptr.Value(10));
      info = double(lusol_obj.get_value(obj.luparm_ptr,10));
    end

    function k = nsing(obj)
      %NSING  number of singularities marked in depcol
      %
      % This method may not be called after updates.  The nsing parameter
      % is only computed after a full factorize.
      %

      % this method only works if no updates have occured.
      obj.update_check();
      %k = double(obj.luparm_ptr.Value(11));
      k = double(lusol_obj.get_value(obj.luparm_ptr,11));
    end

    function k = colrank(obj)
      %RANK  the rank of the matrix determined by the number of independent columns
      %
      % This method uses the LUSOL parameter nsing.  This is only computed
      % after a full factorize.  Thus this method should not be used to
      % determine the rank of a matrix after updates.  In that case look at
      % the flags that are returned by the update methods.
      %
      % Rank determination with LUSOL is more reliable under threshold rook
      % pivoting.  Example options:
      %
      %  options = lusol.luset('pivot','TRP','Ltol1',10)
      %
      n = double(obj.n_ptr.Value);
      nsing = obj.nsing();
      k = n-nsing;
    end

    function d = depcol(obj)
      %DEPCOL  logical vector indicating dependent columns
      %
      % This method may not be called after updates.  The method looks at
      % data that only relevant after a factorize.
      %

      % this method only works if no updates have occured.
      obj.update_check();
      d = obj.depcol_lx;
    end

    function ip = p(obj)
      %P  return row permutation vector
      m = double(obj.m_ptr.Value);
      ip = double(obj.p_ptr.Value(1:m));
    end

    function iq = q(obj)
      %Q  return column permutation vector
      n = double(obj.n_ptr.Value);
      iq = double(obj.q_ptr.Value(1:n));
    end

    % methods to get the matrix factors
    function [L] = L(obj)
       L = [obj.L11; obj.L21]; 
    end
    
    function [U] = U(obj)
        U = [obj.U11, obj.U12t'];
    end
    function [A] = Apq(obj)
        %A = [obj.A11,obj.A12;obj.A21,obj.A22];
        A = obj.A(obj.ap,obj.aq);
    end
    
    function [U11] = getU11(obj)
      %U  get the upper triangular factor U
      %
      % Extract the U factor from LUSOL data and return as a Matlab sparse
      % matrix.
      %
      % Set up:
      %   mylu = lusol(A);
      %
      % Usage:
      %   U1 = mylu.U();
      %   [U2 p q] = mylu.U();
      %   [U2 P Q] = mylu.U('matrix');
      %
      % The first call returns U as a permuted triangle.  The second
      % returns U as an upper triangular matrix with permutation vectors p
      % and q.  The third call returns sparse permutation matrices P and Q.
      % The result of the three calls would produce:
      %   U1(p,q) == P*U1*Q == U2 == upper triangular
      %

      %
      % After a factorize (call to lu1fac) LUSOL stores U by rows at the
      % start of arrays a and indr.  lenr(1:m) stores the number of entries
      % in each row in original order.  locr(1:m) points to the beginning
      % of rows and is stored in original order.
      %
      % Special care must be taken when A is rank deficient.  LUSOL
      % actually stores lenU-nsing entries.  I suppose the extra nsing
      % contained in lenU could be for the zeros on the diagonal.  However,
      % LUSOL seems to handle these implicitly.
      %

      % handle optional function input
      % get basic matrix information
      [m n] = obj.size();
      s = obj.stats();
      
      % obtain required matrix data
      a = obj.a_ptr.Value;
      lenr = obj.lenr_ptr.Value;
      locr = obj.locr_ptr.Value;
      indr = obj.indr_ptr.Value; 
      
      lenU = sum(lenr); 
      % initialize arrays for U triplets
      ui = zeros(lenU,1,'double');
      uj = zeros(lenU,1,'double');
      ua = zeros(lenU,1,'double');

      % array position pointers
      k1 = 1;
      k2 = 1;
      p = obj.p; 
      % loop through (rows?)
      for k = 1:s.nrank1
        % get length of row
        i = p(k); 
        len = double(lenr(i));
        % get location of row
        loc = double(locr(i));
        % set end pointer
        k2 = k1+len-1;
        % load data into triplet arrays
        ui(k1:k2) = i*ones(len,1);
        uj(k1:k2) = double(indr(loc:loc+len-1));
        ua(k1:k2) = double(a(loc:loc+len-1));
        % increment row start pointer
        k1 = k1+len;
      end
      % generate sparse matrix
      U11 = sparse(ui,uj,ua,m,n);
    end
    

    
    function [L0] = L0(obj)
      %L0  get the initial lower triangular factor L0
      %
      % Extracts the initial lower triangular factor from the LUSOL data
      % structure.  LUSOL stores updates to L in product form, thus updates
      % are not included in L0.
      %
      % Set up:
      %   mylu.factorize(A);
      %
      % Usage:
      %   L1 = mylu.get_L0();
      %   [L2 p] = mylu.get_L0();
      %   [L2 P] = mylu.get_L0('matrix');
      %
      % The first call returns L1 as a permuted triangle.  The second and
      % third call will return L2 as a lower triangular matrix with
      % permutation vector p or matrix P.  The result will give:
      %   L1(p,p) == P*L1*P' == L2 == lower triangular
      %

      %
      % After a factorize (call to lu1fac) LUSOL stores non-trivial columns
      % of L at the end of a, indc, and indr.  lenc(1:numL0) stores the
      % number of entries in each column, not including the 1 on the
      % diagonal.  The negatives of the elements of L are stored in a.
      %
      % indc gives the row indices for non-zero elements
      % indr gives the column indices
      %
      % It turns out that off diagonal elements of L are stored in triplet
      % form at the end of a, indc, and indr.  It remains to add the ones
      % on the diagonal.
      %

      % permutation flag, set true if user desires upper triangular U and
      % permutation vectors
      % obtain information from object
      [m n] = obj.size();
      s = obj.stats();
      lena = double(obj.nzmax_ptr.Value);
      % get matrix data
      a = obj.a_ptr.Value;
      indc = obj.indc_ptr.Value;
      indr = obj.indr_ptr.Value;
      % allocate arrays for triplet form
      li = zeros(s.lenL0+m,1);
      lj = zeros(s.lenL0+m,1);
      la = zeros(s.lenL0+m,1);
      % read the triplet form from LUSOL data
      li(1:s.lenL0) = double(indc(lena-s.lenL0+1:lena));
      lj(1:s.lenL0) = double(indr(lena-s.lenL0+1:lena));
      la(1:s.lenL0) = -double(a(lena-s.lenL0+1:lena));
      % add 1's along the diagonal
      li(s.lenL0+1:end) = (1:m)';
      lj(s.lenL0+1:end) = (1:m)';
      la(s.lenL0+1:end) = ones(m,1);
      % create matlab sparse matrix
      L0 = sparse(li,lj,la);
    end
    
    function [inform] = reprow(obj,i1,r)
        nrank = obj.rank;
        r = lusol_obj.vector_process(r,nrank); 
        v_ptr = libpointer('doublePtr', zeros(nrank,1));
        w_ptr = libpointer('doublePtr', zeros(nrank,1));    
        wnew_ptr = libpointer('doublePtr', r);
        
        mode1_ptr = libpointer(obj.int_ptr_class, 1);
        mode2_ptr = libpointer(obj.int_ptr_class, 1);
        irep_ptr = libpointer(obj.int_ptr_class, i1);
        inform_ptr = libpointer(obj.int_ptr_class,0);
        
        calllib('libclusol', 'clu8rpr', ...
            mode1_ptr, ...
            mode2_ptr, ...
            obj.m_ptr, ...
            obj.n_ptr, ...
            irep_ptr, ...
            v_ptr, ...
            w_ptr, ...
            wnew_ptr, ...
            obj.nzmax_ptr, ...
            obj.luparm_ptr, ...
            obj.parmlu_ptr, ...
            obj.a_ptr, ...
            obj.indc_ptr, ...
            obj.indr_ptr, ...
            obj.p_ptr, ...
            obj.q_ptr, ...
            obj.lenc_ptr, ...
            obj.lenr_ptr, ...
            obj.locc_ptr, ...
            obj.locr_ptr, ...
            obj.lenlv_ptr, ...
            obj.li_ptr, ...
            obj.lj_ptr, ...
            obj.lv_ptr, ...
            inform_ptr);
        
        inform = inform_ptr.value; 
        if inform == 7
            error('Update Not Completed- Insufficient Storage. Try increasing nzinit');
        end
    end
    function [inform] = repcol(obj,j1, c)
        nrank = obj.stats.nrank;
        c = lusol_obj.vector_process(c,nrank);
        v_ptr = libpointer('doublePtr', c);
        w_ptr = libpointer('doublePtr', zeros(size(c,1),1));
        
        mode1_ptr = libpointer(obj.int_ptr_class,1);
        mode2_ptr = libpointer(obj.int_ptr_class,1);
        jrep_ptr = libpointer(obj.int_ptr_class,j1);
        inform_ptr = libpointer(obj.int_ptr_class,0);
        diag_ptr = libpointer('doublePtr', 0);
        vnorm_ptr = libpointer('doublePtr', 0); 
        
        calllib('libclusol', 'clu8rpc', ...
            mode1_ptr, ...
            mode2_ptr, ...
            obj.m_ptr, ...
            obj.n_ptr, ...
            jrep_ptr, ...
            v_ptr, ...
            w_ptr, ...
            obj.nzmax_ptr, ...
            obj.luparm_ptr, ...
            obj.parmlu_ptr, ...
            obj.a_ptr, ...
            obj.indc_ptr, ...
            obj.indr_ptr, ...
            obj.p_ptr, ...
            obj.q_ptr, ...
            obj.lenc_ptr, ...
            obj.lenr_ptr, ...
            obj.locc_ptr, ...    
            obj.locr_ptr, ...
            obj.lv_ptr, ...
            obj.li_ptr, ...
            obj.lj_ptr, ...
            obj.lenlv_ptr, ...
            inform_ptr, ...
            diag_ptr, ...
            vnorm_ptr);
        inform = inform_ptr.value;       
    end

    function [inform,d,w1] = swapRows(obj, i1,i2)
        % Swaps rows i1 and nrank+i2 of A and updates the LU factorization
        % accordingly. Returns a rank one update d*w1' that must be applied
        % to U21 to complete the update process. 
  
        nrank = obj.rank; 
        m = size(obj.A,1);
        r1 = obj.A11(i1,:)';
        r2 = obj.A21(i2,:)';
         r12 = obj.A12(i1,:)';
%         %r22 = obj.A22(i2,:)';
         r22 = obj.A(obj.ap(nrank+i2),obj.aq(nrank+1:end))';
         
         r = [r1,r2];
         q = [r12,r22]; 
        
        w1 = lusol_obj.vector_process(r2-r1,nrank);
        v_ptr = libpointer('doublePtr', zeros(nrank,1));
        w_ptr = libpointer('doublePtr', w1);    
        wnew_ptr = libpointer('doublePtr', zeros(nrank,1));
        
        mode1_ptr = libpointer(obj.int_ptr_class, 4);
        mode2_ptr = libpointer(obj.int_ptr_class, 0);
        irep_ptr = libpointer(obj.int_ptr_class, i1);
        inform_ptr = libpointer(obj.int_ptr_class,0);
        
        calllib('libclusol', 'clu8rpr', ...
            mode1_ptr, ...
            mode2_ptr, ...
            obj.m_ptr, ...
            obj.n_ptr, ...
            irep_ptr, ...
            v_ptr, ...
            w_ptr, ...
            wnew_ptr, ...
            obj.nzmax_ptr, ...
            obj.luparm_ptr, ...
            obj.parmlu_ptr, ...
            obj.a_ptr, ...
            obj.indc_ptr, ...
            obj.indr_ptr, ...
            obj.p_ptr, ...
            obj.q_ptr, ...
            obj.lenc_ptr, ...
            obj.lenr_ptr, ...
            obj.locc_ptr, ...
            obj.locr_ptr, ...
            obj.lenlv_ptr, ...
            obj.li_ptr, ...
            obj.lj_ptr, ...
            obj.lv_ptr, ...
            inform_ptr);

        inform = inform_ptr.value;
        if( inform == 7)
            error('lusol:swapRows:lu8rpr', 'Insufficient Storage');
        end        
        
        if (inform == 8)
            % This really shouldn't ever happen
            error('lusol:swapRows:lu8rpr', 'irep out of range'); 
        end
       
        c = sparse(v_ptr.value);
        ei = ((1:m-nrank) == i2)'; 
        %d = sparse(-obj.L21*c - ei);
        
        d = c;     
        for i = obj.lops:1:-1
           d(obj.li(i)) = d(obj.li(i)) - d(obj.lj(i));
        end
        d = obj.B21*(obj.U11 \ d); 
        d = -d - ei;
        % Apply rank one update to U12

        w2 = r22 - r12;
        %obj.U12 = obj.U12 + c*w2';
        %obj.U12t = obj.U12t + w2*c';
        
        % Update S
        %obj.S = obj.S - obj.O*obj.L21*c*w2'; 
         obj.S = obj.S - obj.O*(obj.B21*(obj.U11\ c))*w2';
%         
         lenlv = obj.lenlv_ptr.value; 
         li = double(obj.li_ptr.value);
         lj = double(obj.lj_ptr.value);
         lv = obj.lv_ptr.value;
         lops = sum(lenlv);
         for i=1:lops
           obj.L11(:,lj(i)) = ...
                obj.L11(:,lj(i)) - lv(i)*obj.L11(:,li(i));
%            obj.L21(:,lj(i)) = ...
%                 obj.L21(:,lj(i)) - lv(i)*obj.L21(:,li(i));

           %obj.U12(li(i),:) = ... 
           %         obj.U12(li(i),:) + lv(i)*obj.U12(lj(i),:); 
            %v0 = obj.U12t(:,li(i)) + lv(i)*obj.U12t(:,lj(i)); 
            %obj.U12t(:,li(i)) = v0;
        end
                 
         obj.ap([i1, nrank+i2]) =  obj.ap([nrank+i2, i1]);
         obj.A11(i1,:) = r2;
         obj.A12(i1,:) = r22;
         obj.A21(i2,:) = r1;
         
 
         %obj.A22(i2,:) = r12;
%         
%         %obj.S = obj.S - (obj.O*ei)*w2';
         obj.S = obj.S - (obj.O(:,i2))*w2';
    end
    
    function [inform,s,ej] = swapCols(obj, j1,j2)
        % Swaps columns j1 and nrank+j2 of A and updates the LU factors
        % accordingly 
        
        [~, nrank] = obj.size(); 
        c1 = obj.A11(:,j1);
        c2 = obj.A12(:,j2);
        c = [c1,c2]; 
        c12 = obj.A21(:,j1);
        %c22 = obj.A22(:,j2);
        c22 = obj.A(obj.ap(nrank+1:end),obj.aq(nrank+j2));
%        d = [c12,c22]; 
%        u0 = obj.U12(:,j2);
        %u0 = obj.U12t(j2,:)';
        u0 = obj.solveL(c2); 
        
        %s_old = c22 - obj.L21*u0; 
        s = c22 - obj.A21*obj.solveU(u0);
        % Insert column c1 into position j2
        u = obj.solveL(c1);
        %obj.U12(:,j2) = u;
        %obj.U12t(j2,:) = u';
       % [~,i,v] = find(obj.U12t(j2,:));
        %obj.U12t(j2,i) = 0;
        %[~,i,v] = find(u');
        %obj.U12t(j2,i) = v;
        % Update S
        %obj.S(:,j2) = obj.S(:,j2) - obj.O*obj.L21*(u-u0);
        obj.S(:,j2) = obj.S(:,j2) - obj.O*obj.A21*obj.solveU(u-u0);
        
        v = lusol_obj.vector_process(c2(1:nrank),nrank);
        v_ptr = libpointer('doublePtr', v);
        w_ptr = libpointer('doublePtr', zeros(nrank,1));
        
        mode1_ptr = libpointer(obj.int_ptr_class,1);
        mode2_ptr = libpointer(obj.int_ptr_class,1);
        jrep_ptr = libpointer(obj.int_ptr_class,j1);
        inform_ptr = libpointer(obj.int_ptr_class,0);
        diag_ptr = libpointer('doublePtr', 0);
        vnorm_ptr = libpointer('doublePtr', 0); 
        
        lenLi = obj.stats.lenL;
        calllib('libclusol', 'clu8rpc', ...
            mode1_ptr, ...
            mode2_ptr, ...
            obj.m_ptr, ...
            obj.n_ptr, ...
            jrep_ptr, ...
            v_ptr, ...
            w_ptr, ...
            obj.nzmax_ptr, ...
            obj.luparm_ptr, ...
            obj.parmlu_ptr, ...
            obj.a_ptr, ...
            obj.indc_ptr, ...
            obj.indr_ptr, ...
            obj.p_ptr, ...
            obj.q_ptr, ...
            obj.lenc_ptr, ...
            obj.lenr_ptr, ...
            obj.locc_ptr, ...    
            obj.locr_ptr, ...
            obj.lv_ptr, ...
            obj.li_ptr, ...
            obj.lj_ptr, ...
            obj.lenlv_ptr, ...
            inform_ptr, ...
            diag_ptr, ...
            vnorm_ptr);
        inform = inform_ptr.value;
        
        if( inform == 7)
            error('lusol:swapCols', 'Insufficient Storage');
        end        
       
        
        lenlv = obj.lenlv_ptr.value; 
        li = double(obj.li_ptr.value);
        lj = double(obj.lj_ptr.value);
        lv = obj.lv_ptr.value;
        lops = sum(lenlv(1:2));
        
        obj.lops = lops;
        obj.li = li;
        obj.lj = lj;
        obj.lv = lv; 
        
        lenLf = obj.stats.lenL;
        if lenLf - lenLi ~= lops
            error('lusol:swapCols', 'Invalid L update');
        end
        
        for i=1:lops
            obj.L11(:,lj(i)) = ...
                obj.L11(:,lj(i)) - lv(i)*obj.L11(:,li(i));
            
%             obj.L21(:,lj(i)) = ...
%                 obj.L21(:,lj(i)) - lv(i)*obj.L21(:,li(i));
%              
            %obj.U12(li(i),:) = ... 
            %        obj.U12(li(i),:) + lv(i)*obj.U12(lj(i),:); 
           % obj.U12t(:,li(i)) = obj.U12t(:,li(i)) + lv(i)*obj.U12t(:,lj(i));
        end
 
         obj.aq([j1 nrank+j2]) = obj.aq([nrank+j2 j1]);        
         obj.S(:,j2) = obj.S(:,j2) + obj.O*(c12-c22);
         obj.A11(:,j1) = c2;
         obj.A12(:,j2) = c1;
         obj.A21(:,j1) = c22;
%         obj.A22(:,j2) = c12;
         ej = ((1:nrank) == j1)';
    end
    
    

    
    function [inform] =  swapFac(obj, a_r, a_c, s_r, s_c)
        nrank = obj.rank;
        [m,n] = size(obj.A);
        inform = 0; 
        obj.U11 = obj.getU11(); 
        obj.B21 = obj.A21;
        

        if a_c <= nrank
            [inform,s,ej] = obj.swapCols(a_c, s_c);
        end
 
        if a_r <= nrank
            [inform,d,w1] = obj.swapRows(a_r, s_r);
        end
           
        if(obj.stats.nrank ~= nrank || inform == 2)
           error('lusol:swapFac', 'Singular U11');
        end
        
%         ej1 = ((1:nrank) == a_c)';
%         ej2 = ((1:n-nrank) == s_c)'; 
%         ei1 = ((1:nrank) == a_r)';
%         ei2 = ((1:m-nrank) == s_r)'; 
%         
%         E11 = [c(:,2)-c(:,1), ei1]; 
%         F11 = [ej1'; r(:,2)'-r(:,1)']; 
%         
%         E22 = [d(:,1)-d(:,2),ei2]; 
%         F22 = [ej2'; q(:,1)'-q(:,2)'];
%         
%         E12 = [c(:,1)-c(:,2),ei1];
%         F12 = [ej2' ; q(:,2)'-q(:,1)']; 
%         
%         E21 = [d(:,2)-d(:,1), ei2]; 
%         F21 = [ej1'; r(:,1)'-r(:,2)'];
%         
%         
% 
%         D12 = obj.solveA(E12);
%         D12 = D12*F12;
%         D21 = obj.solveAt(F21')*E21';
%         M = obj.solveA(E11)*(eye(2) - F11*obj.solveA(E11))^-1*(obj.solveAt(F11')');
%       
% 
%         f1 = obj.O*E22*F22;
%         f2 = obj.O*A21*M*A12;
%         f3 = obj.O*D21'*A12; 
%         f4 = obj.O*A21*D12 ;
%         f5   = obj.O*E21*F21*D12;
% 
%         obj.S = obj.S +f1+ f2 - f3 - f4- f5;
       if a_c <= nrank
           l = sparse(obj.solveUt(ej));
            %obj.L21 = obj.L21 + s*l'; 
            %obj.S = obj.S - obj.O*s*l'*obj.U12;
            
            %obj.S = obj.S - obj.O*s*(obj.U12t*l)';
            obj.S = obj.S -obj.O*s*(obj.solveLt(l)'*obj.A12); 
       end
        
       if a_r <= nrank
            l2 = sparse(obj.solveUt(w1)); 
            %obj.L21 = obj.L21 + d*l2';
            %obj.S = obj.S - obj.O*d*l2'*obj.U12;
            
            %obj.S = obj.S  -obj.O*d*(obj.U12t*l2)';           
            obj.S = obj.S -obj.O*d*(obj.solveLt(l2)'*obj.A12); 
       end


    end
    function [err] = updateA12(obj,a_r, a_c, s_r, s_c)
        nrank = obj.stats.nrank; 
        m = size(obj.A,1); 

        new_row = obj.A(obj.ap(nrank+s_r), obj.aq(nrank+1:end));
        old_row = obj.A12(a_r,:);
        %new_row(s_c) = obj.A(obj.ap(nrank+s_r),obj.aq(a_c)); 
        w2 = new_row - old_row; 
        ei1 = ((1:nrank)' == a_r);
        ei2 = ((1:m-nrank)' == s_r); 
        %c1  = obj.solveL(ei1);
        %c2 = - obj.A21*obj.solveU(c1) -ei2; 
        c2 = -obj.A21*obj.solveA(ei1) - ei2; 
        obj.S = obj.S + (obj.O*c2)*w2;
        obj.A12(a_r,:) = new_row; 
%        obj.A22(s_r,:) = obj.A22(s_r,:) - w2;

        % Update Col s_c
        new_col = obj.A(obj.ap(1:nrank),obj.aq(a_c));
        new_col(a_r) = obj.A(obj.ap(nrank+s_r),obj.aq(a_c));
        v1 = new_col- obj.A12(:,s_c);
        c1 = obj.solveL(v1); 
        c2 = - obj.A21*obj.solveU(c1);
        obj.S(:,s_c) = obj.S(:,s_c) + obj.O*c2; 
        obj.A12(:,s_c) = new_col;
        
%         S = obj.O*(obj.A22 - obj.A21*obj.A11^-1*obj.A12);
%         err = norm(S-obj.S,'fro')/norm(S,'fro');
        err = 0;
    end
    
    function [err] = updateA21(obj,a_r,a_c,s_r,s_c)
        nrank = obj.stats.nrank;
        n = size(obj.A, 2); 

        %Update row s_r
        new_row = obj.A(obj.ap(a_r),obj.aq(1:nrank));
        w1 = new_row - obj.A21(s_r,:); 
        %l = obj.solveUt(w1'); 
        %c1 = obj.solveLt(l); 
        c1 = obj.solveAt(w1');
        obj.S = obj.S - obj.O(:,s_r)* (c1'*obj.A12); 
        obj.A21(s_r,:) = new_row; 

         
        new_col = obj.A(obj.ap(nrank+1:end),obj.aq(nrank+s_c));
        new_col(s_r) = obj.A(obj.ap(a_r),obj.aq(nrank+s_c));
        old_col = obj.A21(:,a_c);
        v2 = new_col - old_col; 
        ej1 = (1:nrank)'== a_c; 
        ej2 = (1:n-nrank) == s_c; 
%         l = obj.solveUt(ej1); 
%         obj.S = obj.S - (obj.O*v2)*(obj.solveLt(l)'*obj.A12 ) ;        
        c2 = obj.solveAt(ej1);
        obj.S = obj.S - (obj.O*v2)*(c2'*obj.A12);
        obj.S(:,s_c) = obj.S(:,s_c) - obj.O*v2; 
        obj.A21(:,a_c) = new_col; 

        %obj.A22(:,s_c) = obj.A22(:,s_c) - v2;
       
%         S = obj.O*(obj.A22 - obj.A21*obj.A11^-1*obj.A12);
%         err = norm(S-obj.S,'fro')/norm(S,'fro');
        err = 0;
    end
    
    function updateA22(obj, a_r, a_c, s_r, s_c)
        
        nrank = obj.stats.nrank; 
   
           % Update row s_r
        %new_row= obj.A(obj.ap(a_r),obj.aq(nrank+1:end));
   
        %new_row(s_c) = obj.A(obj.ap(a_r),obj.aq(a_c)); 
        %old_row = obj.A(obj.ap(nrank+s_r),obj.aq(nrank+1:end));
        %old_row(s_c) = obj.A(obj.ap(nrank+s_r),obj.aq(a_c)); 
        %w2 = new_row - old_row; 
        %obj.S = obj.S+ obj.O(:,s_r)*w2; 
        
        % Update col s_c
%         new_col = obj.A(obj.ap(nrank+1:end),obj.aq(a_c)); 
%         new_col(s_r) = obj.A(obj.ap(a_r),obj.aq(a_c)); 
%         old_col = obj.A(obj.ap(nrank+1:end), obj.aq(nrank+s_c)); 
%         old_col(s_r) = obj.A(obj.ap(a_r), obj.aq(nrank+s_c)); 
%         v2 = new_col - old_col; 
%         obj.S(:,s_c) = obj.S(:,s_c) + obj.O*v2; 
        %obj.A22(:,s_c) = new_col; 
        

       % obj.A22(s_r,:) = new_row; 
        
        
    end
    
    function  [err] = updateA11(obj, a_r,a_c,s_r,s_c)
        nrank = obj.stats.nrank; 
        new_col = obj.A(obj.ap(1:nrank),obj.aq(nrank+s_c)); 
        old_col = obj.A11(:,a_c);
        obj.A11(:, a_c) = new_col; 
        
        new_row = obj.A(obj.ap(nrank+s_r),obj.aq(1:nrank));
        new_row(a_c) = obj.A(obj.ap(nrank+s_r),obj.aq(nrank+s_c)); 
        old_row = obj.A11(a_r,:);   
        obj.A11(a_r,:) = new_row; 
   

        v1 = (new_col - old_col);
        u1 = new_row - old_row; 
        ei = ((1:nrank)' == a_r); 

        ej = ((1:nrank) == a_c); 

        
        cj = obj.solveAt(ej');
        x1 = obj.solveAt(u1');
        obj.repcol(a_c,new_col);
        obj.reprow(a_r,new_row);
        y1 = obj.solveA(v1);
        ci = obj.solveA(ei); 
          %E11 = [new_col - old_col , ei];       
         %F11 = [ej; new_row-old_row];
%          vnorm = vecnorm(E11(:,1));
%          unorm = vecnorm(F11(2,:));
%          E11(:,1) = E11(:,1)/sqrt(vnorm);
%          E11(:,2) = E11(:,2)*sqrt(unorm);
%         F11(2,:) = F11(2,:)/sqrt(unorm);
%          F11(1,:) = F11(1,:)*sqrt(vnorm);
%          F11(2,:) = F11(2,:)/sqrt(unorm);
%          F11(1,:) = F11(1,:)*sqrt(vnorm);
%         N1 = obj.solveA(E11);
%         N2 = obj.solveAt(F11')';
%         C = eye(2) - F11*obj.solveA(E11);

%         M2 = (N1*(C\N2)); 
%         obj.S = obj.S + obj.O*obj.A21*M2*obj.A12;    
        y2 = obj.A21*y1;
        ci2 = obj.A21*ci;
        cj2 = cj'*obj.A12;
        x2 = x1'*obj.A12;
        obj.S = obj.S + (obj.O*y2)*cj2 + (obj.O*ci2)*x2;
%         S = obj.O*(obj.A22 - obj.A21*obj.A11^-1*obj.A12);
%         err = norm(S-obj.S,'fro')/norm(S,'fro');
        err = 0;
    end
    
    function swapFac2(obj, a_r,a_c,s_r,s_c)
        nrank = obj.stats.nrank;
        
        if(a_r ~= nrank+1 && a_c == nrank+1)
            obj.swapFacRows(a_r,s_r);
  
        
        elseif (a_r == nrank+1 && a_c ~= nrank+1)
            obj.swapFacCols(a_c,s_c)
                
        elseif (a_r == nrank+1 && a_c == nrank+1)
            err = 0;
            return
        else
            err1 = obj.updateA12(a_r,a_c,s_r,s_c);
            err2 = obj.updateA21(a_r,a_c,s_r,s_c);
            %obj.updateA22(a_r,a_c,s_r,s_c);
            err3 = obj.updateA11(a_r,a_c,s_r,s_c);
        
            obj.ap([a_r,nrank+s_r]) = obj.ap([nrank+s_r,a_r]);
            obj.aq([a_c,nrank+s_c]) = obj.aq([nrank+s_c,a_c]);
        end

    end
    
    function swapFac3(obj, a_r,a_c,s_r,s_c, alpha,beta)
         [m,n] = size(obj.A); 
         nrank = obj.stats.nrank;
        if(a_r ~= nrank+1 && a_c == nrank+1)
            obj.swapFacRows(a_r,s_r);
  
        
        elseif (a_r == nrank+1 && a_c ~= nrank+1)
            obj.swapFacCols(a_c,s_c)
                
        elseif (a_r == nrank+1 && a_c == nrank+1)
            return
        else
            
            a2 = obj.A12(:,s_c);
            a4 = obj.A(obj.ap(nrank+1:end),obj.aq(nrank+s_c));
            
            b3= obj.A21(s_r,:)'; 
            b4 = obj.A(obj.ap(nrank+s_r),obj.aq(nrank+1:end))';
            ej1 = (1:nrank)' ==a_c;
            ej2 = (nrank+1:n)' == nrank+s_c;        
            ei1 = (1:nrank)' == a_r;
            ei2 = (nrank+1:m)' == nrank+s_r;
            
            s11 = -alpha; 
            s22 = ej1'*obj.solveA(ei1);
            w = obj.solveAt(b3); 
            s12 = w(a_r);
            
            v = obj.solveA(a2); 
            s21 = v(a_c);
            s = [s11, s12; s21, s22]; 
            v1 = obj.A21*v -a4; 
            v2 = obj.A21*obj.solveA(ei1) + ei2; 
            v2_ = v2 + s12/alpha*v1;
            v2_(s_r) = 1; 
            
            
            w1 = w'*obj.A12 - b4'; 
            w2  = obj.solveAt(ej1)'*obj.A12 + ej2'; 
            w2_ = w2 + s21/alpha*w1;
            w2_(s_c) = 1; 
            s_ = [alpha, 0 ; 0, beta]; 
            %E = ((obj.O*[v1,v2])/s)*[w1; w2];
            E = obj.O*[v1,v2_]*[-1/alpha*w1 ; 1/beta*w2_];
            obj.S = obj.S + E; 
            
            %Update A11
            A11 = obj.A11;
            new_col = obj.A(obj.ap(1:nrank),obj.aq(nrank+s_c)); 
            old_col = obj.A11(:,a_c);
            obj.A11(:, a_c) = new_col; 
        
            new_row = obj.A(obj.ap(nrank+s_r),obj.aq(1:nrank));
            new_row(a_c) = obj.A(obj.ap(nrank+s_r),obj.aq(nrank+s_c)); 
            old_row = obj.A11(a_r,:);   
            obj.A11(a_r,:) = new_row; 


            v1 = (new_col - old_col);
            u1 = new_row - old_row; 
            ei = ((1:nrank)' == a_r); 
            ej = ((1:nrank) == a_c); 
            
            obj.repcol(a_c,new_col);
            obj.reprow(a_r,new_row);
            
            % Update perm and A12/A21
            obj.ap([a_r, nrank+s_r]) = obj.ap([nrank+s_r,a_r]);
            obj.aq([a_c, nrank+s_c]) = obj.aq([nrank+s_c,a_c]);
            
            obj.A12(a_r,:) = obj.A(obj.ap(a_r),obj.aq(nrank+1:end));
            obj.A12(:,s_c) = obj.A(obj.ap(1:nrank),obj.aq(nrank+s_c));
            
            obj.A21(:,a_c) = obj.A(obj.ap(nrank+1:end),obj.aq(a_c)); 
            obj.A21(s_r,:) = obj.A(obj.ap(nrank+s_r),obj.aq(1:nrank)); 
            
           end
    end
    function swapFacRows(obj,a_r,s_r)
        nrank = obj.rank; 
        m = size(obj.A,1);
        %Update A12
        new_row = obj.A(obj.ap(nrank+s_r),obj.aq(nrank+1:end)); 
        w2 = new_row - obj.A12(a_r,:);
        ei1 = ((1:nrank)' == a_r);
        ei2 = ((1:m-nrank)' == s_r); 
        c2 = -obj.A21*obj.solveA(ei1) - ei2; 
        obj.S = obj.S + obj.O*c2*w2;
        obj.A12(a_r,:) = new_row; 
%         obj.A22(s_r,:) = obj.A22(s_r,:) - w2;
        
        %Update A21
        new_row = obj.A(obj.ap(a_r),obj.aq(1:nrank));
        w1 = new_row - obj.A21(s_r,:);        
        c1 = obj.solveAt(w1');
        obj.S = obj.S - obj.O(:,s_r)* (c1'*obj.A12); 
        obj.A21(s_r,:) = new_row;   
        
        %Update A11
        new_row = obj.A(obj.ap(nrank+s_r),obj.aq(1:nrank));
        w1 = new_row - obj.A11(a_r,:);
        ei = ((1:nrank)' == a_r);
        n1 = obj.solveA(ei);
        n2 = obj.solveAt(w1'); 
        c = 1 + w1*n1;
        obj.S = obj.S + (obj.O*(obj.A21*n1))*c^-1*(n2'*obj.A12);
        obj.reprow(a_r,new_row);
        obj.A11(a_r,:) = new_row;
        obj.ap([a_r,nrank+s_r]) = obj.ap([nrank+s_r,a_r]);
        
    end
    
    function swapFacCols(obj,a_c,s_c)
        nrank = obj.rank;
        
        %Update A12 
        new_col = obj.A(obj.ap(1:nrank),obj.aq(a_c));
        v1 = new_col- obj.A12(:,s_c);
        c2 = obj.solveA(v1);
        obj.S(:,s_c) = obj.S(:,s_c) - obj.O*obj.A21*c2; 
        obj.A12(:,s_c) = new_col;     
        
        %Update A21
        new_col = obj.A(obj.ap(nrank+1:end),obj.aq(nrank+s_c));
        v2 = new_col - obj.A21(:,a_c); 
        ej1 = (1:nrank)'== a_c;     
        c2 = obj.solveAt(ej1);
        obj.S = obj.S - (obj.O*v2)*(c2'*obj.A12);
        obj.S(:,s_c) = obj.S(:,s_c) - obj.O*v2; 
        obj.A21(:,a_c) = new_col; 
%         obj.A22(:,s_c) = obj.A22(:,s_c) - v2;
        
        
        %Update A11
        new_col = obj.A(obj.ap(1:nrank),obj.aq(nrank+s_c));
        v1 = new_col - obj.A11(:,a_c);
        ej = (1:nrank)' == a_c;
        n1 = obj.solveA(v1);
        n2 = obj.solveAt(ej);
        c = 1 + ej'*n1;
        obj.S = obj.S + (obj.O*(obj.A21*n1))*c^-1*(n2'*obj.A12);
        obj.repcol(a_c,new_col);
        obj.A11(:,a_c) = new_col;
        obj.aq([a_c,nrank+s_c]) = obj.aq([nrank+s_c,a_c]);       
        
    end
    function [alpha, s_r, s_c] = maxS(obj)

        [~, s_c] = max(sqrt(sum(obj.S.^2)));
        nrank = obj.stats.nrank;
        %max_col = obj.A22(:,s_c) - obj.L21*obj.U12(:,s_c);
        %max_col_old = obj.A(obj.ap(nrank+1:end),obj.aq(nrank+s_c)) - obj.L21*(obj.U12t(s_c,:))';
        %max_col = obj.A(obj.ap(nrank+1:end),obj.aq(nrank+s_c)) - obj.A21*obj.solveU(obj.U12t(s_c,:)');
        max_col = obj.A(obj.ap(nrank+1:end),obj.aq(nrank+s_c)) - ...
                            obj.A21*obj.solveU(obj.solveL(obj.A12(:,s_c)));
            
        [~, s_r] = max(abs(max_col));
        alpha = full(max_col(s_r));
        
  
    end

    
    function [beta, a_r, a_c] = maxA11inv(obj,alpha, s_r,s_c)
        %% Estimates the maximum entry of A11^-1 where A11 is the matrix 
        % formed by rows [1:nrank,nrank+s_r] and cols [1:nrank,nrank+s_c]
        % of A
       
        %% Compute the factorization A11 = L11*U11 
        % L11 = [L  , 0]    U11 = [U,   u12]
        %       [l21', 1]          [0, alpha]
        % L and U are the nrank x nrank submatrix of the 
        % current truncated LU factorization

        nrank = obj.stats.nrank; 
        %l21_old = obj.L21(s_r,:)';
        l21 = obj.solveUt(obj.A21(s_r,:)');
        %u12 = obj.U12(:,s_c);
        %u12 = obj.U12t(s_c,:)'; 
        u12 = obj.solveL(obj.A12(:,s_c));

        %% Compute B =  Omega*A11^-1
        % [v1,v2] = [Omega1,Omega2]*U11^-1
        Omega1 = randn(20,nrank); 
        Omega2 = randn(20,1); 
        v1 = obj.solveUt(Omega1');
        v2 = 1/alpha*(Omega2' - u12'*v1)'; 
        
        %Compute B = ([v1',v2])/L11;
        B2 = v2;
        B1 = obj.solveLt(v1-l21*v2'); 
        B = [B1', B2];
        
        %% Find maximum element
        [~,a_r] = max(sqrt(sum(B.^2)));
        e_mcol = zeros(nrank+1,1);
        e_mcol(a_r,1) = 1;
       
        %Compute u = L11\max_col;
        u1 = obj.solveL(e_mcol(1:nrank)); 
        u2 = e_mcol(nrank+1) - l21'*u1; 
        u = [u1 ; u2]; 
        v2 = u(nrank+1)/alpha;
        v1 = obj.solveU(u(1:nrank) - v2*u12); 
        max_col = [v1 ; v2];
        [~,a_c] = max(abs(max_col)); 
        beta = full(max_col(a_c));

    end
    
    function [nswaps,fval] = srlu(obj, f, maxswaps)
        nswaps = 0;
        errors = zeros(1,maxswaps);
        is = zeros(1,maxswaps);
        js = zeros(1,maxswaps); 
        nrank = obj.rank;
        for i = 1:maxswaps
             if(mod(i,100) == 0)
                obj.refactor('pivot', 'TCP', 'nzinit', 5*10^6);
            end
            [alpha, s_r, s_c] = obj.maxS(); 
            if abs(alpha) < 10^-15
                break
            end
      
            [beta, a_r, a_c] = obj.maxA11inv(alpha,s_r, s_c);

            fi = abs(alpha*beta);
            fval= fi; 
            if fi < f
                break            
            end 
            udiag = obj.diagU(); 
            det1 = sum(log(abs(udiag)));
            
            obj.swapFac3(a_r,a_c,s_r,s_c,alpha,beta);
            

            nswaps = nswaps+1;
            
            udiag = obj.diagU();
            det2 = sum(log(abs(udiag))); 
            if(det2-det1 < log(f)-.001)
                error('lusol:srlu', 'det not increasing enough');
            end
        end
       %Compute the submatrices
       nrank = obj.rank;
       obj.L11 = obj.mulL(eye(nrank));
       Linv  = obj.solveL(eye(nrank));
       obj.U11 = obj.getU11();
       obj.L21 = obj.A21 / obj.U11; 
       obj.U12 = Linv*obj.A12; 
    end
    
    function [M] = cur(obj)
       %A = [obj.A11, obj.A12 ; obj.A21 , obj.A22];
       A = obj.Apq();
       L = [obj.L11 ; obj.L21];
       U = [obj.U11, obj.U12t']; 
       
       nrank = obj.stats.nrank;
       d = obj.U11(sub2ind([nrank,nrank],obj.p,obj.q));
       D = diag(d);       
       Udt = U' ./ d; 
       M1 = (Udt'*Udt) \ (Udt'*A');
       M1 = (D\ M1)';
       
       M = (L'*L)\ (L'*M1);
       
    end
    
    function [condL, condU] = cond(obj)
       L = [obj.L11;obj.L21];
       U = [obj.U11 , obj.U12t']; 
       nrank = obj.stats.nrank;
       
       d = obj.U11(sub2ind([nrank,nrank],obj.p,obj.q));    
       Udt = U' ./ d;       
       [~,Rl]  = qr(L,0);
       [~,Ru] = qr(Udt,0);
       
       condL = cond(full(Rl));
       condU = cond(full(Ru));
       
    end
            
    function [U1,S0,V1] = qrcur(obj)
        nrank = obj.stats.nrank;
        L = [obj.L11; obj.L21]; U = [obj.getU11(), obj.U12t']; 
        m = obj.m;
        n = obj.n;
        apinv = zeros(m,1);aqinv = zeros(n,1); apinv(obj.ap) = 1:m; aqinv(obj.aq) = 1:n; 
        L = L(apinv,:); U = U(:,aqinv);      
        
        [~,R1] = qr(L,0);
        Ql = L/R1;
        C1 = Ql'*obj.A;
        
        [~,R2] = qr(U',0);
        Qu = U'/R2;
        C2 = C1*Qu; 
        
        [U0,S0,V0] = svds(C2,nrank);
        U1 = Ql*U0;
        V1 = Qu*V0;
        
    end
    
    
    function [C,M,R] = stable_cur(obj, krank)
        nrank = obj.stats.nrank;
        C = obj.A(obj.ap(1:end),obj.aq(1:nrank));
        R = obj.A(obj.ap(1:nrank),obj.aq(1:end)); 
        
        [Qc,R1] = qr(C,0);
        C1 = Qc'*obj.Apq; 
      %  [C1,R1 ] = qr(C, obj.Apq, 0);

        
        %[Qr,R2] =qr(R',0);
        [C2, R2] = qr(R',C1',0); 

        
        %B = Qc'*obj.A(obj.ap,obj.aq)*Qr; 
        B = C2';
        obj.B = B; 
        %B = R1'\(C'*obj.A(obj.ap,obj.aq)*R')/R2; 
        [U0,S0,V0] = svds(B,krank);
        Bk = U0*S0*V0'; 
        
        M = R1\(Bk) / R2'; 
      % obj.Qc = Qc;
        obj.Bk = Bk;
      %  obj.Qr = Qr'; 
        
        %obj.Rc = R1;
       % obj.Rr = R2; 
    end

    function [diagU] = diagU(obj)
       [m,n] = obj.size;
       nrank = obj.stats.nrank;
       diagU = zeros(m,1);
       diagU_ptr = libpointer('doublePtr', diagU);
       ret_inform_ptr = libpointer(obj.int_ptr_class,0);
       calllib('libclusol', 'clu9diagu',...
          obj.m_ptr,...
          obj.n_ptr,...
          obj.nzmax_ptr,...
          obj.luparm_ptr,...
          diagU_ptr,...
          obj.a_ptr,...
          obj.indc_ptr,...
          obj.indr_ptr,...
          obj.p_ptr,...
          obj.q_ptr,...
          obj.lenr_ptr,...
          obj.locr_ptr,...
          ret_inform_ptr);

       diagU = diagU_ptr.value;
       diagU = diagU(1:nrank);
    end 
    
    
    function [err] = L2error(obj)
        L = [obj.L11; obj.L21];
        U = [obj.getU11(), obj.U12];
        %U = [obj.getU11(), obj.U12];
        %Apq = [obj.A11,obj.A12; obj.A21, obj.A22];
        Apq = obj.Apq();
        %M = (L \ Apq) / U;
        %E = Apq - L*M*U;
        E = Apq - L*U;        
        err = norm(E,'fro')/norm(obj.A,'fro');
    end
    
    function [err] = cur_error(obj)
        nrank = obj.stats.nrank;
        %C = obj.A(obj.ap(1:end),obj.aq(1:nrank));
       % R = obj.A(obj.ap(1:nrank),obj.aq(1:end)); 
        
        %P1 = obj.Qc'*obj.Apq*obj.Qr'; 
        %P = obj.Rc'\(C'*obj.Apq*R')/obj.Rr; 
        err = sqrt(1-norm(obj.B,'fro')^2/norm(obj.A,'fro')^2 + norm(obj.B-obj.Bk,'fro')^2/norm(obj.A,'fro')^2);
    end
    

    function [angles] = angle_error(obj)
        L = [obj.L11 ; obj.L21];
        [Ul,~,~] = svds(L, 5);
        Apq = obj.Apq(); 
        [Ua,~,~]= svds(Apq, 5); 
        angles = svd(Ul'*Ua); 
    end
    function [err] = L2errorB(obj,block_size)
        L = [obj.L11; obj.L21];
        %U = [obj.getU11(), obj.U12t'];
        U = [obj.getU11(), obj.U12t'];
        %Apq = [obj.A11,obj.A12; obj.A21, obj.A22];
        Apq = obj.Apq();
        %M = (L \ Apq) / U;
        %E = Apq - L*M*U;
        err = 0; 
        n = size(Apq,2);
        blocks = [0:block_size:n,n];
        for i=1:length(blocks)-1
            E = Apq(:,blocks(i)+1:blocks(i+1)) - L*U(:,blocks(i)+1:blocks(i+1));
            
            err = err+ norm(E,'fro')^2;
        end
        err = err/norm(obj.A,'fro')^2;
        err = sqrt(err); 
        
    end

    function [err] = CURerrorB(obj,block_size)
        L = [obj.L11; obj.L21];
        U = [obj.getU11(), obj.U12t'];

        Apq = obj.Apq();
        M = obj.cur(); 
        err = 0; 
        n = size(Apq,2);
        blocks = [0:block_size:n,n];
        for i=1:length(blocks)-1
            E = Apq(:,blocks(i)+1:blocks(i+1)) - L*M*U(:,blocks(i)+1:blocks(i+1));
            
            err = err+ norm(E,'fro')^2;
        end
        err = err/norm(obj.A,'fro')^2;
        err = sqrt(err); 
        
    end
   
    function [err] = facerror(obj)
        err = [obj.facerror11(),obj.facerror12(),obj.facerror21()];
        err = full(err);
    end
    function [err] = facerror11(obj)
        nrank = obj.rank;
        U11 = obj.getU11();
        L11 = obj.mulL(eye(nrank));
        err = max(max(abs(obj.A11 - L11*U11)));
    end
    
    function [err] = facerror12(obj)
        nrank = obj.stats.nrank;
        [~,n] = size(obj.A);
        err = max(max(abs(obj.A12 - obj.L11*obj.U12t')));
        %err = 1;
    end
 
    function [err] = facerror21(obj)
        n = obj.n; 
        U11  = obj.getU11();
        err = max(max(abs(obj.A21 - obj.L21*U11)));
        %err = 1; 
    end
    
    function [err] = serr(obj)
        nrank = obj.stats.nrank;
        obj.A22 = obj.A(obj.ap(nrank+1:end),obj.aq(nrank+1:end));
        S0 = obj.O*obj.A22 - obj.O*obj.A21*obj.A11^-1*obj.A12; 
        err = max(max(abs(S0-obj.S)));
    end
    function [e,s,U0] = svd_err(obj,sk,use_m)
        [U0,~,~] = svds(obj.A,sk); 
        [m,n] = size(obj.A); 
        L = [obj.L11; obj.L21];
        U = [obj.U11, obj.U12t'];
        
        apinv = zeros(m,1);
        aqinv = zeros(n,1); 
        apinv(obj.ap) = 1:m; 
        aqinv(obj.aq) = 1:n; 
        
        L = L(apinv,:); 
        U = U(:,aqinv); 
        if use_m
            [U1, S, ~] = svds(L*U, sk);
        else
            M = obj.cur();
             [U1, S, ~] = svds(L*M*U, sk);
        end
        e = svd(U0'*U1); 
        s = diag(S); 
    end
    % solve methods

    function [X inform resid] = solve(obj,B,mode)
      %SOLVE  solve systems with matrix factors
      %
      % This function solves all of the relavent systems of equations.  If
      % right hand side B is a matrix, it will solve for matrix X.
      %
      % Input:
      %  B = right hand size, can be a vector or a matrix
      %  mode = solution mode, see table below
      %
      % Output:
      %  X = solution matrix
      %  inform = status flag vector, one element for each column of B
      %  resid = 1-norm of residuals for each solve
      %
      % Modes:
      %  1    X  solves   L * X = B
      %  2    X  solves   L'* X = B
      %  3    X  solves   U * X = B
      %  4    X  solves   U'* X = B
      %  5    X  solves   A * X = B (default)
      %  6    X  solves   A'* X = B
      %
      % inform flags:
      %  0 = successful solve
      %  1 = if U is singular, and residual is non-zero
      %

      % set default mode
      if nargin < 3
        mode = 5;
      end

      % get size of factorized matrix
      [m n] = obj.size();

      % get size of B
      [Br Bc] = size(B);

      % compute size of X
      Xc = Bc;
      switch mode
        case 1 % X  solves   L X = B
          % B must have m rows
          if Br ~= m, error('lusol:solve','B has incorrect size.'); end
          % X is m by Bc
          Xr = m;
        case 2 % X  solves   L'X = B
          % B must have m rows
          if Br ~= m, error('lusol:solve','B has incorrect size.'); end
          % X is m by Bc
          Xr = m;
        case 3 % X  solves   U X = B
          % B must have m rows
          if Br ~= m, error('lusol:solve','B has incorrect size.'); end
          % X is n by Bc
          Xr = n;
        case 4 % X  solves   U'X = B
          % B must have n rows
          if Br ~= n, error('lusol:solve','B has incorrect size.'); end
          % X is m by Bc
          Xr = m;
        case 5 % X  solves   A X = B
          % B must have m rows
          if Br ~= m, error('lusol:solve','B has incorrect size.'); end
          % X is n by Bc
          Xr = n;
        case 6 % X  solves   A'X = B
          % B must have n rows
          if Br ~= n, error('lusol:solve','B has incorrect size.'); end
          % X is m by Bc
          Xr = m;
        otherwise
          error('lusol:solve','unrecognized mode.')
      end

      % allocate space for output X
      X = zeros(Xr,Xc);
      inform = zeros(1,Bc);
      resid = zeros(1,Bc);

      % compute solutions for all columns
      for j = 1:Bc
        [X(:,j) inform(j) resid(j)] = obj.clu6sol(B(:,j),mode);
      end

    end

    function [X inform resid] = solveA(obj,B)
      %SOLVEA  solve A*X = B.
      %
      % See also: lusol.solve
      [X inform resid] = obj.solve(B,5);
    end

    function [X inform resid] = solveAt(obj,B)
      %SOLVEAT  solve A'*X = B.
      %
      % See also: lusol.solve
      [X inform resid] = obj.solve(B,6);
    end

    function [X inform] = solveL(obj,B)
      %SOLVEL  solve L*X = B.
      %
      % See also: lusol.solve
      [X inform] = obj.solve(B,1);
    end
   
    function [X inform] = solveL11(obj, B)
      nrank = obj.stats.nrank;
      [m,~] = obj.size;
      Br = size(B,1);
      if Br ~= nrank
          error('lusol:solveL11', 'B has incorrect size');
      end
      B = [B ; zeros(m-nrank,size(B,2))];
      [X inform] = obj.solve(B,1);
      X = X(1:nrank,:); 
    end
    
    function [X inform] = solveLt(obj,B)
      %SOLVELT  solve L'*X = B.
      %
      % See also: lusol.solve
      [X inform] = obj.solve(B,2);
    end
    
    function [X inform] = solveL11t(obj,B)
       nrank = obj.stats.nrank;
       [m,~] = obj.size;
       Br = size(B,1);
       if Br ~= nrank
           error('lusol:solveL11', 'B has incorrect size');
       end       
       B = [B ; zeros(m-nrank,size(B,2))];
       [X inform] = obj.solve(B,2);
       X = X(1:nrank,:); 
    end

    function [X inform resid] = solveU(obj,B)
      %SOLVEU  solve U*X = B.
      %
      % See also: lusol.solve
      [X inform resid] = obj.solve(B,3);
    end
    
    function [X inform] = solveU11(obj, B)
      nrank = obj.stats.nrank; 
      [m,~] = obj.size(); 
      B = [B ; zeros(m-nrank, size(B,2))];  
      
      [X inform resid] = obj.solve(B,3);
      X = X(1:nrank,:); 
    end

    function [X inform resid] = solveUt(obj,B)
      %SOLVEUT  solve U'*X = B.
      %
      % See also: lusol.solve
      [X inform resid] = obj.solve(B,4);
    end

    function [X inform] = solveU11t(obj,B)
      nrank = obj.stats.nrank; 
      [~,n] = obj.size(); 
      B = [B ; zeros(n-nrank, size(B,2))]; 
      [X inform ~] = obj.solve(B,4); 
      X = X(1:nrank,:); 
        
    end
    % multiply methods

    function Y = mul(obj,X,mode)
      %MUL  compute matrix multiplies with various factors
      %
      % Perform matrix-vector or matrix-matrix multiply with desired matrix
      % factors.
      %
      % Usage:
      %  Y = mylu.mul(X,mode)
      %
      % Input:
      %  X = matrix or vector to compute multiply
      %  mode = multiply mode, see table below
      %
      % Output:
      %  Y = product of desired factor and X
      %
      % Mode:
      %  1    Y = L * X
      %  2    Y = L'* X
      %  3    Y = U * X
      %  4    Y = U'* X
      %  5    Y = A * X (default)
      %  6    Y = A'* X
      %

      % set default mode
      if nargin < 3
        mode = 5;
      end

      % get size of X
      [Xr Xc] = size(X);

      % get size of factored matrix
      [m n] = obj.size();

      % compute size of X
      Yc = Xc;
      switch mode
        case 1 %  Y = L X
          % X must have m rows
          if Xr ~= m, error('lusol:mul','X has incorrect size.'); end
          % Y is m by Xc
          Yr = m;
        case 2 %  Y = L'X
          % X must have m rows
          if Xr ~= m, error('lusol:mul','X has incorrect size.'); end
          % Y is m by Xc
          Yr = m;
        case 3 %  Y = U X
          % X must have n rows
          if Xr ~= n, error('lusol:mul','X has incorrect size.'); end
          % Y is m by Xc
          Yr = m;
        case 4 %  Y = U'X
          % X must have m rows
          if Xr ~= m, error('lusol:mul','X has incorrect size.'); end
          % Y is n by Xc
          Yr = n;
        case 5 %  Y = A X
          % X must have n rows
          if Xr ~= n, error('lusol:mul','X has incorrect size.'); end
          % Y is m by Xc
          Yr = m;
        case 6 %  Y = A'X
          % X must have m rows
          if Xr ~= m, error('lusol:mul','X has incorrect size.'); end
          % Y is n by Xc
          Yr = n;
        otherwise
          error('lusol:mul','unrecognized mode.')
      end

      % allocate space for output Y
      Y = zeros(Yr,Yc);

      % compute solutions for all columns
      for j = 1:Xc
        [Y(:,j)] = obj.clu6mul(X(:,j),mode);
      end

    end

    function Y = mulA(obj,X)
      %MULA  compute Y = A*X.
      %
      % See also: lusol.mul
      Y = obj.mul(X,5);
    end
    
    function Y = mulA22(obj, X)
       nrank = obj.stats.nrank;
       X = [zeros(nrank,size(X,2)) ; X];
       Y = obj.mul(X,5);
       Y = Y(nrank+1:end,:); 
    end

    function Y = mulAt(obj,X)
      %MULAT  compute Y = A'*X.
      %
      % Warning: this does not seem to work at the moment.
      %
      % See also: lusol.mul
      Y = obj.mul(X,6);
    end
    
    function Y = mulA22t(obj, X)
       nrank = obj.stats.nrank;
       X = [zeros(nrank,size(X,2)) ; X]; 
       Y =obj.mul(X,6);
       Y = Y(nrank+1:end,:); 
    end

    function Y = mulL(obj,X)
      %MULL  compute Y = L*X.
      %
      % See also: lusol.mul
      Y = obj.mul(X,1);
    end
    
    function Y = mulL11(obj,X)
        nrank = obj.stats.nrank;
        [m, ~]= obj.size; 
        Xr = size(X,1);
        if Xr ~= nrank
            error("lusol:mulL11", "X has incorrect size");
        end
        
        X = [X ; zeros(m-nrank, size(X,2))];
        Y = obj.mul(X,1);
        Y = Y(1:nrank,:);
    end
    
    function Y = mulL21(obj, X)
      nrank = obj.stats.nrank;
      [m, ~] = obj.size;
      Xr = size(X,1);
      if Xr ~= nrank
          error("lusol:mulL21", "X has incorrect size");
      end
      
      X = [X ; zeros(m-nrank,size(X,2))];
      Y = obj.mul(X,1);
      Y = Y(nrank+1:end,:); 
    end
    
    function Y = mulLt(obj,X)
      %MULLT  compute Y = L'*X.
      %
      % See also: lusol.mul
      Y = obj.mul(X,2);
    end
    
    function Y = mulL21t(obj, X)
        nrank = obj.stats.nrank;
        [m,~] = obj.size; 
        Xr = size(X,1);
        if Xr ~= m-nrank
            error('lusol:mulL21',"X has incorrect size");
        end
        X = [zeros(nrank,size(X,2)); X];
        Y = obj.mul(X,2);
        Y = Y(1:nrank,:); 
    end
    function Y = mulU(obj,X)
      %MULU  compute Y = U*X.
      %
      % See also: lusol.mul
      Y = obj.mul(X,3);
    end
    
    function Y = mulU11(obj, X)
       [~,n] = obj.size; 
       nrank = obj.stats.nrank;
       Xr = size(X,1);
       if Xr ~= nrank
           error('lusol:mulU11', "X has incorrect size");
       end
       X = [X ; zeros(n-nrank,size(X,2))];
       Y = obj.mul(X,3);
       Y = Y(1:nrank,:);        
    end
    
    function Y = mulU12(obj, X)
       [~,n] = obj.size; 
       nrank = obj.stats.nrank;
       Xr = size(X,1);
       if Xr ~= n-nrank
           error('lusol:mulU12', "X has incorrect size");
       end
       X = [zeros(nrank,size(X,2)); X];
       Y = obj.mul(X,3);
       Y = Y(1:nrank,:); 
    end
    
    function Y = mulUt(obj,X)
      %MULUT  compute Y = U'*X.
      %
      % See also: lusol.mul
      Y = obj.mul(X,4);
    end

    function Y = mulU12t(obj, X)
       nrank = obj.stats.nrank; 
       [m,n] = obj.size; 
       X = [X ; zeros(m-nrank,size(X,2))]; 
       Y = obj.mul(X,4);
       Y = Y(nrank+1:end,:);
    end
  

    function [inform,v] = r1mod(obj,v,w,beta)
      %R1MOD  update LU factorization to perform rank-1 update A+beta*v*w'
      %
      % Usage:
      %  inform = mylu.r1mod(v,w,beta)
      %
      % Input:
      %  v = vector of length m
      %  w = vector of length n
      %  beta = scalar value
      %
      % Outputs:
      %  inform = status flag
      %
      % On exit:
      %  inform =  -1 if the rank of U decreased by 1.
      %  inform =  0  if the rank of U stayed the same.
      %  inform =  1  if the rank of U increased by 1.
      %  inform =  7  if the update was not completed (lack of storage).
      %

      locr = obj.locr_ptr.Value;
      lenr = obj.lenr_ptr.Value;
      indr = obj.indr_ptr.Value;
      lv1 = locr(265);
      lenv = lenr(265);
      lv2 = lv1 + lenv - 1;
      jvs = indr(lv1:lv2);
      
      
      
      % get size of matrix
      [m n] = obj.size();

      % handle variable input
      if nargin < 4 || isempty(beta)
        beta = 1.0;
      end
      
      % check and process input vectors
      v = lusol_obj.vector_process(v,m);
      w = lusol_obj.vector_process(w,n);

      % prepare temporary libpointers
      beta_ptr = libpointer('doublePtr',beta);
      v_ptr = libpointer('doublePtr',v);
      w_ptr = libpointer('doublePtr',w);
      ret_inform_ptr = libpointer(obj.int_ptr_class,0);

      % call library function
      calllib('libclusol','clu9mod', ...
        obj.m_ptr, ...
        obj.n_ptr, ...
        beta_ptr, ...
        v_ptr, ...
        w_ptr, ...
        obj.nzmax_ptr, ...
        obj.luparm_ptr, ...
        obj.parmlu_ptr, ...
        obj.a_ptr, ...
        obj.indc_ptr, ...
        obj.indr_ptr, ...
        obj.p_ptr, ...
        obj.q_ptr, ...
        obj.lenc_ptr, ...
        obj.lenr_ptr, ...
        obj.locc_ptr, ...
        obj.locr_ptr, ...
        ret_inform_ptr);

      
      % prepare function output
      inform = ret_inform_ptr.Value;
      v = v_ptr.Value; 
    end
    
    function [inform,lmax] = luclear(obj, c)
      [m,~] = obj.size();
      lmax_ptr = libpointer('doublePtr', 0); 
      ret_inform_ptr = libpointer(obj.int_ptr_class,0);
      c = lusol_obj.vector_process(c, m); 
      c_ptr = libpointer('doublePtr', c); 
      calllib('libclusol', 'clu9clr', ...
        obj.m_ptr, ...
        obj.n_ptr, ...
        obj.nzmax_ptr, ...
        obj.luparm_ptr, ...
        obj.parmlu_ptr, ...
        obj.a_ptr, ...
        obj.indc_ptr, ...
        obj.indr_ptr, ...
        obj.p_ptr, ...
        obj.q_ptr, ...
        obj.lenr_ptr, ...
        obj.locc_ptr, ...
        obj.locr_ptr, ...
        c_ptr, ...
        lmax_ptr, ...
        ret_inform_ptr); 

      inform = ret_inform_ptr.Value; 
      lmax = lmax_ptr.Value;
    end
    % inherited, not implemented methods

    function [L, U, p,q] = lupq(obj)
        L = obj.L;
        U = obj.U;
        p = obj.ap;
        q = obj.aq;
    end
    function TF = eq(obj,x)
      %EQ  Not implemented, returns false
      TF = false;
    end

    function TF = ge(obj,x)
      %GE  Not implemented, returns false
      TF = false;
    end

    function TF = gt(obj,x)
      %GT  Not implemented, returns false
      TF = false;
    end

    function TF = le(obj,x)
      %LE  Not implemented, returns false
      TF = false;
    end

    function TF = lt(obj,x)
      %LT  Not implemented, returns false
      TF = false;
    end

    function TF = ne(obj,x)
      %NE  Not implemented, returns false
      TF = false;
    end

    function addlistener(obj,varargin)
      %ADDLISTENER  Not implemented
    end

    function notify(obj,varargin)
      %NOTIFY  Not implemented
    end

  end

end % classdef
