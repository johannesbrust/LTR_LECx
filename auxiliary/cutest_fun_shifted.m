function [varargout] = cutest_fun_shifted( x, varargin )
  % Evaluate objective function
  % Usage:       f = obj(x)      evaluates function value only
  %          [f,g] = obj(x)  evaluates function value and gradient
  %        [f,g,H] = obj(x)  evaluates function value, gradient and Hessian

  % Global "shift" parameter
  global rho;
  
  if nargout > 3
      error( 'obj: too many output arguments' );
  end

  if nargin == 1
      % Compute objective function value
      if nargout > 1 
         % Gradient is requested
        [varargout{1}, varargout{2}] = cutest_obj(x);
        varargout{1} = varargout{1} + rho/2*(x'*x);
        varargout{2} = varargout{2} + rho*x;
      else
        varargout{1} =  cutest_obj(x) + rho/2*(x'*x);
      end
  elseif nargin == 2  
      % Only gradient is requested
      [~, varargout{1}] = cutest_obj(x);
      
      varargout{1} = varargout{1} + rho*x;
      
      % Using "cutest_grad" may cause crashes for some problems
      %[varargout{1}] = cutest_grad(x);
  end