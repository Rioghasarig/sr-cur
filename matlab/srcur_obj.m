classdef srcur_obj < handle
    properties (Access=private)
        corelu = 0; 
        rank = 0; 

        O = 0; 
        m = 0;
        n = 0;
        
        ap = [];
        aq = []; 
        
        S = [];
        A = []; 
        A11 = [];
        A12 = [];
        A21 = [];
        A22 = []; 
        B = [];
        Bk = []; 
        nzinit2 = 0; 
    end
    
    methods (Static)
        function options = srcur_options(varargin)
              % get the input parser
              in_parse = inputParser;

              % storage parameters
              in_parse.addParameter('nzinit',0,@(x) x>=0);
              in_parse.addParameter('nzinit2',0,@(x) x>=0);
              in_parse.addParameter('pivot','TCP',@(x) ismember(x,{'TRP','TCP'}));
              in_parse.addParameter('rank',1, @(x) x>=0);
              
              in_parse.parse(varargin{:}); 
              options = in_parse.Results;
        end
        
        
        
    end
    
   
    methods (Access = public)
        function obj = srcur_obj(A, varargin)
            obj.A = A; 
            obj.m = size(A,1);
            obj.n = size(A,2); 
            options = obj.srcur_options(varargin{:}); 
            
            
            trlu = lusol(A, 'pivot', options.pivot, 'rank', options.rank, 'nzinit', options.nzinit );
            nrank = trlu.stats.nrank; 
            obj.rank = nrank; 

            obj.ap = trlu.p;
            obj.aq = trlu.q; 
            
            L = trlu.getL0();
            L = L(obj.ap,obj.ap(1:nrank)); 
            U = trlu.getU();
            U = U(obj.ap(1:nrank),obj.aq); 
            L21 = L(nrank+1:end,1:nrank);       
            U12 = U(1:nrank,nrank+1:end);
      
            obj.O = randn(20,obj.m-nrank);
            OLU = obj.O*L21;
            OLU = OLU*U12;
            obj.S = obj.O*obj.A(obj.ap(nrank+1:end),obj.aq(nrank+1:end)) - ...
                OLU;
            
            
            obj.A11 = obj.A(obj.ap(1:nrank),obj.aq(1:nrank));
            obj.A12 = obj.A(obj.ap(1:nrank),obj.aq(nrank+1:end));
            obj.A21 = obj.A(obj.ap(nrank+1:end),obj.aq(1:nrank));
            obj.A22 = obj.A(obj.ap(nrank+1:end),obj.aq(nrank+1:end));
            
            % Refactor A11
            obj.nzinit2 = options.nzinit2; 
            obj.corelu = lusol(obj.A11, 'pivot', 'TCP', 'rank', nrank, 'nzinit', options.nzinit2); 
            
      
        end
        function refactor(obj)

          nrank = obj.rank;
          obj.corelu = lusol(obj.A11, 'pivot', 'TCP', 'rank', nrank, 'nzinit', obj.nzinit2); 
        
          if( obj.corelu.stats.nrank ~= nrank)
              error('lusol:refactor','A11 not full rank');
          end
          
          L11 = obj.corelu.getL0();
          U11 = obj.corelu.getU(); 
          
          U12 = (L11 \ obj.A12);
          L21 = obj.A21 / U11; 

          OLU = obj.O*L21;
          OLU = OLU*U12;
          obj.S = obj.O*obj.A(obj.ap(nrank+1:end),obj.aq(nrank+1:end)) - ...
                    OLU;

        end

        function [alpha, s_r, s_c] = maxS(obj)

            [~, s_c] = max(sqrt(sum(obj.S.^2)));
            nrank = obj.rank;
            max_col = obj.A(obj.ap(nrank+1:end),obj.aq(nrank+s_c)) - ...
                                obj.A21*obj.corelu.solveU(obj.corelu.solveL(obj.A12(:,s_c)));

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

            nrank = obj.rank; 
            %l21_old = obj.L21(s_r,:)';
            l21 = obj.corelu.solveUt(obj.A21(s_r,:)');
            %u12 = obj.U12(:,s_c);
            %u12 = obj.U12t(s_c,:)'; 
            u12 = obj.corelu.solveL(obj.A12(:,s_c));

            %% Compute B =  Omega*A11^-1
            % [v1,v2] = [Omega1,Omega2]*U11^-1
            Omega1 = randn(20,nrank); 
            Omega2 = randn(20,1); 
            v1 = obj.corelu.solveUt(Omega1');
            v2 = 1/alpha*(Omega2' - u12'*v1)'; 

            %Compute B = ([v1',v2])/L11;
            B2 = v2;
            B1 = obj.corelu.solveLt(v1-l21*v2'); 
            B = [B1', B2];

            %% Find maximum element
            [~,a_r] = max(sqrt(sum(B.^2)));
            e_mcol = zeros(nrank+1,1);
            e_mcol(a_r,1) = 1;

            %Compute u = L11\max_col;
            u1 = obj.corelu.solveL(e_mcol(1:nrank)); 
            u2 = e_mcol(nrank+1) - l21'*u1; 
            u = [u1 ; u2]; 
            v2 = u(nrank+1)/alpha;
            v1 = obj.corelu.solveU(u(1:nrank) - v2*u12); 
            max_col = [v1 ; v2];
            [~,a_c] = max(abs(max_col)); 
            beta = full(max_col(a_c));

        end
        
        function [beta, a_r,a_c] =  maxA11tinv(obj, alpha,s_r,s_c)
            nrank = obj.rank;
            Omega = rand(20,nrank+1);
            
            a21 = obj.A(obj.ap(nrank+s_r),obj.aq(1:nrank))'; 
            a12 = obj.A(obj.ap(1:nrank),obj.aq(nrank+s_c));
            B11 = obj.A(obj.ap([1:nrank,nrank+s_r]),obj.aq([1:nrank,nrank+s_c]));
            A11 = obj.A(obj.ap([1:nrank]),obj.aq(1:nrank));
            
            w21 = obj.corelu.solveAt(a21);
            w12 = obj.corelu.solveA(a12);
            
            X3 = [eye(nrank),w12; zeros(1,nrank),1];
            X2 = [obj.A(obj.ap([1:nrank]),obj.aq([1:nrank])), zeros(nrank,1);zeros(1,nrank),alpha];
            X1 = [eye(nrank),zeros(nrank,1); w21',1]; 
            
            Y1 = [eye(nrank),-w21;zeros(1,nrank),1];
            Y3 = [eye(nrank),zeros(nrank,1);-w12',1];
            
            O = Omega;
            Omega(:,nrank+1) = Omega(:,nrank+1) - Omega(:,1:nrank)*w21;
            Omega(:,1:nrank) = obj.corelu.solveA(Omega(:,1:nrank)')';
            Omega(:,nrank+1) = Omega(:,nrank+1)/alpha; 
            Omega(:,1:nrank) = Omega(:,1:nrank) - Omega(:,nrank+1)*w12'; 
            
            [~,a_c] =  max(sum(Omega.^2)); 
            e = zeros(nrank+1,1);
            e(a_c) = 1; 
            
            e(nrank+1) = e(nrank+1)- w12'*e(1:nrank);
            e(1:nrank) = obj.corelu.solveAt(e(1:nrank)); 
            e(nrank+1) = e(nrank+1)/alpha;
            e(1:nrank) = e(1:nrank) - w21*e(nrank+1);
            
            [~,a_r] = max(abs(e)); 
            beta = full(e(a_r)); 
            
        end
        function swapFacRows(obj,a_r,s_r)
            nrank = obj.rank; 
            m = size(obj.A,1);
            %Update A12
            new_row = obj.A(obj.ap(nrank+s_r),obj.aq(nrank+1:end)); 
            w2 = new_row - obj.A12(a_r,:);
            ei1 = ((1:nrank)' == a_r);
            ei2 = ((1:m-nrank)' == s_r); 
            c2 = -obj.A21*obj.corelu.solveA(ei1) - ei2; 
            obj.S = obj.S + obj.O*c2*w2;
            obj.A12(a_r,:) = new_row; 
    %         obj.A22(s_r,:) = obj.A22(s_r,:) - w2;

            %Update A21
            new_row = obj.A(obj.ap(a_r),obj.aq(1:nrank));
            w1 = new_row - obj.A21(s_r,:);        
            c1 = obj.corelu.solveAt(w1');
            obj.S = obj.S - obj.O(:,s_r)* (c1'*obj.A12); 
            obj.A21(s_r,:) = new_row;   

            %Update A11
            new_row = obj.A(obj.ap(nrank+s_r),obj.aq(1:nrank));
            w1 = new_row - obj.A11(a_r,:);
            ei = ((1:nrank)' == a_r);
            n1 = obj.corelu.solveA(ei);
            n2 = obj.corelu.solveAt(w1'); 
            c = 1 + w1*n1;
            obj.S = obj.S + (obj.O*(obj.A21*n1))*c^-1*(n2'*obj.A12);
            obj.corelu.reprow(a_r,new_row);
            obj.A11(a_r,:) = new_row;
            obj.ap([a_r,nrank+s_r]) = obj.ap([nrank+s_r,a_r]);

        end

        function swapFacCols(obj,a_c,s_c)
            nrank = obj.rank;

            %Update A12 
            new_col = obj.A(obj.ap(1:nrank),obj.aq(a_c));
            v1 = new_col- obj.A12(:,s_c);
            c2 = obj.corelu.solveA(v1);
            obj.S(:,s_c) = obj.S(:,s_c) - obj.O*obj.A21*c2; 
            obj.A12(:,s_c) = new_col;     

            %Update A21
            new_col = obj.A(obj.ap(nrank+1:end),obj.aq(nrank+s_c));
            v2 = new_col - obj.A21(:,a_c); 
            ej1 = (1:nrank)'== a_c;     
            c2 = obj.corelu.solveAt(ej1);
            obj.S = obj.S - (obj.O*v2)*(c2'*obj.A12);
            obj.S(:,s_c) = obj.S(:,s_c) - obj.O*v2; 
            obj.A21(:,a_c) = new_col; 
    %         obj.A22(:,s_c) = obj.A22(:,s_c) - v2;


            %Update A11
            new_col = obj.A(obj.ap(1:nrank),obj.aq(nrank+s_c));
            v1 = new_col - obj.A11(:,a_c);
            ej = (1:nrank)' == a_c;
            n1 = obj.corelu.solveA(v1);
            n2 = obj.corelu.solveAt(ej);
            c = 1 + ej'*n1;
            obj.S = obj.S + (obj.O*(obj.A21*n1))*c^-1*(n2'*obj.A12);
            obj.corelu.repcol(a_c,new_col);
            obj.A11(:,a_c) = new_col;
            obj.aq([a_c,nrank+s_c]) = obj.aq([nrank+s_c,a_c]);       

        end

        function swapFac(obj, a_r,a_c,s_r,s_c, alpha,beta)
             [m,n] = size(obj.A); 
             nrank = obj.rank;
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
            s22 = ej1'*obj.corelu.solveA(ei1);
            w = obj.corelu.solveAt(b3); 
            s12 = w(a_r);
            
            v = obj.corelu.solveA(a2); 
            s21 = v(a_c);
            s = [s11, s12; s21, s22]; 
            v1 = obj.A21*v -a4; 
            v2 = obj.A21*obj.corelu.solveA(ei1) + ei2; 
            v2_ = v2 + s12/alpha*v1;
            v2_(s_r) = 1; 
            
            
            w1 = w'*obj.A12 - b4'; 
            w2  = obj.corelu.solveAt(ej1)'*obj.A12 + ej2'; 
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
            
            obj.corelu.repcol(a_c,new_col);
            obj.corelu.reprow(a_r,new_row);
            
            % Update perm and A12/A21
            obj.ap([a_r, nrank+s_r]) = obj.ap([nrank+s_r,a_r]);
            obj.aq([a_c, nrank+s_c]) = obj.aq([nrank+s_c,a_c]);
            
            obj.A12(a_r,:) = obj.A(obj.ap(a_r),obj.aq(nrank+1:end));
            obj.A12(:,s_c) = obj.A(obj.ap(1:nrank),obj.aq(nrank+s_c));
            
            obj.A21(:,a_c) = obj.A(obj.ap(nrank+1:end),obj.aq(a_c)); 
            obj.A21(s_r,:) = obj.A(obj.ap(nrank+s_r),obj.aq(1:nrank)); 
            
           end
        end
        
        function [nswaps,fval] = srlu(obj, f, maxswaps)
            nswaps = 0;
            for i = 1:maxswaps
                 if(mod(i,50) == 0)
                    obj.refactor();
                 end
                
                [alpha, s_r, s_c] = obj.maxS(); 
                if abs(alpha) < 10^-15
                    break
                end

                [beta, a_r, a_c] = obj.maxA11inv(alpha,s_r, s_c);
                [beta2,a_r2,a_c2] = obj.maxA11tinv(alpha,s_r,s_c);
                fi = abs(alpha*beta);
                fval= fi; 
                if fi < f
                    break            
                end 
                udiag = obj.corelu.diagU(); 
                det1 = sum(log(abs(udiag)));

                obj.swapFac(a_r,a_c,s_r,s_c,alpha,beta);


                nswaps = nswaps+1;

                udiag = obj.corelu.diagU();
                det2 = sum(log(abs(udiag))); 
                if(det2-det1 < log(f)-.001)
                    error('lusol:srlu', 'det not increasing enough');
                end
            end
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
        function [err] = cur_error(obj)
            nrank = obj.stats.nrank;
            %C = obj.A(obj.ap(1:end),obj.aq(1:nrank));
            % R = obj.A(obj.ap(1:nrank),obj.aq(1:end)); 

            %P1 = obj.Qc'*obj.Apq*obj.Qr'; 
            %P = obj.Rc'\(C'*obj.Apq*R')/obj.Rr; 
            err = sqrt(1-norm(obj.B,'fro')^2/norm(obj.A,'fro')^2 + norm(obj.B-obj.Bk,'fro')^2/norm(obj.A,'fro')^2);
        end

        function [err] = facerror11(obj)
            nrank = obj.rank;
            U11 = obj.corelu.getU();
            L11 = obj.corelu.mulL(eye(nrank));
            err = max(max(abs(obj.A11 - L11*U11)));
        end
        
        function [err] = serr(obj)
            nrank = obj.rank;
            obj.A22 = obj.A(obj.ap(nrank+1:end),obj.aq(nrank+1:end));
            S0 = obj.O*obj.A22 - obj.O*obj.A21*obj.A11^-1*obj.A12; 
            err = max(max(abs(S0-obj.S)));
        end
    end
end
    
    