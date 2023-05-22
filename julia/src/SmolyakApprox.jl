
# This module provides a Fortran translation of the codes
# written by Kenneth L. Judd, Lilia Maliar, Serguei Maliar 
# and Rafael Valero for their paper (2014),  "Smolyak method 
# for solving dynamic economic models: Lagrange interpolation,  
# anisotropic grid and adaptive domain" Journal of Economic Dynamics  
# and Control 44, 92 123 (henceforth, JMMV (2014))

# Copyright 2014 by Lilia Maliar, Serguei Maliar and Rafael Valero. All 
# rights reserved. The code may be used, modified and redistributed under  
# the terms provided in the file "License_Agreement.txt".

# translation to Fortran: Moritz Lenel
# translation to Julia from Fortran: Joseph Kupferberg 

# Smolyak_Elem_Isotrop is a routine that constructs the subindices of the
# Smolyak elements (grid points and basis functions) for the isotropic case;
# see "Smolyak method for solving dynamic economic models: Lagrange interpo-
# lation, anisotropic grid and adaptive domain" by Kenneth L. Judd, Lilia Maliar, 
# Serguei Maliar and Rafael Valero, (2014), Journal of Economic Dynamics and 
# Control 44, 92 123 (henceforth, JMMV (2014)), Section 2.2.3
#
# This version: Novenber 5, 2014. First version: December 17, 2012.
# This translation: February 2020
# -------------------------------------------------------------------------
# Input:   "d" is the number of dimensions (the number of state  variables) 
#          "mu" is the level of approximation
#  
# Output:  "Smolyak_elem_iso" is the vector of subindices of unidimensional 
#           elements (Smolyak grid points or polynomial basis function) that 
#           constitute a multidimensional element (Smolyak grid point or 
#           polynomial basis function) for the isotropic case

function Smolyak_Elem_Isotrop(d, mu) 
                            
                            

    # 1. Identify the indices of disjoint sets A's, i1,...,id, that satisfy the 
    # Smolyak rule, d<=i1+i2+...+id<=|i|; see equation (1) in JMMV(2014)
    # -------------------------------------------------------------------------
       
        Smol_rule = [];     # This will be a matrix of subindices i1,...,id of 
                            # disjoint sets A_i's such that their sum across all  
                            # dimensions, |i|=i1+i2+...+id,  satisfies  the 
                            # Smolyak rule, d<=|i|<=d+mu; initially, the matrix 
                            # is empty; at the intermediate steps j=0,...,mu, 
                            # the subindices satisfy d<=|i|=d+j
                           
        incr_Smol_rule = Int.(ones(1,d));
                            # The matrix of subindices for j=0 is a 
                            # 1-by-d vector of ones, (1,...,1)
         
        
                            # The matrix of subindices i1,...,id of disjoint   
                            # sets A_i's, such that their sum across all dimensions,  
                            # |i|=i1+i2+...+id, is exactly equal to d+j, i.e., 
                            # |i|=d+j; initially, the matrix is empty; when j 
                            # increases, this matrix is concatinated to matrix  
                            # "Smol_rule" obtained for the previous value of j,
                            # i.e., j-1
       for j = 0:mu  
           
         prev_incr = copy(incr_Smol_rule); 
                           # Ffor the previous value of j, call "incr_Smol_rule"  
                           # as "prev_incr"
         
         # Identify new subindices of unidimensional sets A's i1,i2,...,id that 
         # jointly satisfy the Smolyak rule as the sum of subindices increases 
         # from j-1 to j; see JMMV (2014), Section 2.2 
         #---------------------------------------------------------------------                
         if (j > 0)
            
           m = size(prev_incr,1);
                                       # Compute the number of columns in the
                                       # previous matrix of subindices that
                                       # satisfy the Smolyak rule "prev_incr"
           incr_Smol_rule = [];  # Initially, the matrix of subindices is 
                                            # empty   
           aux = Int.(zeros(m,d));           # Allocate memory to an initial auxiliary
                                       # matrix that will be added to
                                       # "prev_incr" 
            for id = 1:d
                aux_new = copy(aux);         # New auxiliary matrix is equal to the old 
                                       # one
                aux_new[:,id] .= 1;      # For a dimension i, set elements of this
                                       # new auxiliary matrix to 1 instead of 0
                                       
                augmented = prev_incr .+ aux_new; 
                                       # Increase the subinices of
                                       # "prevoius_incr" by 1 in dimension id
                if isempty(incr_Smol_rule)
                    incr_Smol_rule = copy(augmented)
                else
                    incr_Smol_rule = [incr_Smol_rule; augmented];
                                       # Concatenate "incr_Smol_rule" and
                                       # "augmented": the first row of "augmented"
                                       # goes after the last row of
                                       # "incr_Smol_rule" 
                end
            end
           
            
         end
        # TODO: turned off sorting for now
         # incr_Smol_rule = sortslices(unique(incr_Smol_rule,dims = 1),dims=1)
                                       # Eliminate the repeated indices in the
                                       # matrix "incr_Smol_rule"   
        incr_Smol_rule = unique(incr_Smol_rule,dims = 1)
         # Concatenate the matrix of newly constructed indices to the previously 
         # obtained matrix of indices satisfying the Smolyak rule
         #---------------------------------------------------------------------
        
         if isempty(Smol_rule)
            Smol_rule = copy(incr_Smol_rule); # E.g., for mu=1 and d=2,
         else
            Smol_rule = [Smol_rule; incr_Smol_rule]           # Smol_rule=[1 1; 1 2; 2 1]
         end
       end
      
       n_comb = size(Smol_rule,1);
                  # The total number of combinations of indices, i1,...,id, 
                  # of the unidimensional disjoint sets A's that satisfy the 
                  # Smolyak rule; e.g., for mu=1 and d=2, n_comb=3 such as (1)
                  # i1=1 and i2=1, (2) i1=1 and i2=2, (3) i1=2 and i2=1
                  
    
    # 2. Construct the multidimensional indices of elements as a Cartesian product
    # of unidimensional indices
    # -------------------------------------------------------------------------
    
    Smolyak_elem_iso = []; 
                    # The matrix of multidimensional indices of elements (points)
                    #  belonging to the disjoint subsets A's that satisfy the 
                    # Smolyak rule (every single point is indexed); initially, 
                    # the matrix is empty
                  
    for i = 1:n_comb            # For each combination of subindices i1,...,id 
                                # that satisfies the Smolyak rule
        
       
        incr_indices = [];      # This is a matrix of multidimensional indices
                                # of unidimensional grid points that will be
                                # added to the final multidimensional matrix of 
                                # indices "Smolyak_elem_iso", as n_comb increases;
                                # initially, the matrix is empty
    
        one_comb = Smol_rule[i,:]; 
                                # Consider the i-th row of "Smol_rule" (i.e.,
                                # one particular combination of the subindices 
                                # of the disjoint sets satisfying the Smolyak
                                # rule); e.g., for mu=1 and d=2, Smol_rule=
                                # [1 1; 1 2; 2 1] and one_comb for i=1 is [1 1];
                                # 1-by-d
        for jd = 1:d            # For each dimension jd, ...
            
            prev_indices = incr_indices;   
                
            # Compute the indices of elements (points) of the unidimensional 
            # set
            #----------------------------------------------------------------
            if one_comb[jd] == 1                        
                               # Take a jd-th element of the row vector "one_comb" 
                               # that corresponds to dimension jd (this is the 
                               # subindex of the unidimensional disjoint set A_i 
                               # from which this element comes from; if an element
                               # (point) is from the disjoint set
                               # A_1,...              
                indices_elem_jd = 1;                        
                               # A_1 contains one element; this element is indexed "1"
       
            elseif one_comb[jd] == 2                    
                               # If an element (point) is from the disjoint set
                               # A_2,...
                indices_elem_jd = [2,2^(one_comb[jd]-1)+1];            
                               # A_2 contains two elements; these elements are
                               # indexed "2" and "3", so that indices_elem_jd=[2;3]
            else
                indices_elem_jd = (2^(one_comb[jd]-2)+2 : 2^(one_comb[jd]-1)+1); 
                               # The subsequent disjoint sets contain the elements 
                               # from m(one_comb[jd]-1)+1=2^(one_comb[jd]-2)+2 to 
                               # m(one_comb[jd])=2^(one_comb[jd]-1)+1; e.g., the 
                               # disjoint set A_4 contains the elements indexed 
                               # 6, 7, 8 and 9, so that indices_elem_jd=[6,7,8,9]
            end                      
                               
               # Create a Cartesian product of two sets, "prev_indices" and 
               # "indices_elem_jd"
               #-----------------------------------------------------------
               
            a = copy(prev_indices);   # Call one set "a" 
            b = copy(indices_elem_jd);# Call the other set "b"
            z = [];         # Initially, the Cartesian product is empty
            if (isempty(b))     # If "b" is empty, ...
                z = a;          # The Cartesian product "z" is given by "a"
            elseif (isempty(a)) # If "a" is empty, ...
                z = b;          # The Cartesian product "z" is given by "b"   
            else                # For the case of non-empty sets, ...
                
                a_columns = size(a,1);  
                b_columns = size(b,1);
                if b_columns == 1
                    z = hcat(a,Int.(ones(a_columns,1)*b))
                else
                    for k = 1:b_columns;
                        if isempty(z)
                            z = hcat(a,Int.(ones(a_columns,1))*b[k,:]); 

                        else    
                            z = vcat(z,hcat(a, Int.(ones(a_columns,1))*b[k,:]));
                        end 
                    end
                end
            end
            incr_indices = copy(z);   #                
        end

        if isempty(Smolyak_elem_iso)
            Smolyak_elem_iso = copy(incr_indices); 
                            # Construct the matrix of multidimensional indices 
                            # of elements by concatenating "incr_indices" to 
                            # the previous "Smolyak_elem_iso" obtained for the 
                            # previous combination "n_comb"
        else
            Smolyak_elem_iso = [Smolyak_elem_iso;incr_indices]; 
        end
    end
    return(Smolyak_elem_iso)
end 


# Smolyak_Elem_Anisotrop is a routine that selects a subset of the subindices 
# of Smolyak elements corresponding to the given anisotropic case from a set of 
# subindices of the Smolyak isotropic elements; see "Smolyak method for solving  
# dynamic economic models: Lagrange interpolation, anisotropic grid and   
# adaptive domain" by Kenneth L. Judd, Lilia Maliar, Serguei Maliar and  
# Rafael Valero, (2014), Journal of Economic Dynamics and Control 44, 92�123 
# (henceforth, JMMV (2014)), Section 2.2.3 
#
# This version: Novenber 5, 2014. First version: December 17, 2012.
# This translation: February 2020
# -------------------------------------------------------------------------
# Input:   "smol_elem_iso"         is the matrix of subindices of unidi- 
#                                  mensional elements that constitute multi-
#                                  dimensional elements for the isotropic case
#          "vector_mus_dimensions" is the vector of the levels of
#                                  approximations in all dimensions  
# Output:  "smol_elem_ani"         is the matrix of subindices of unidimen-
#                                  sional elements that constitute multidi-
#                                  mensional elements for the given anisotropic 
#                                  case
# 
# -------------------------------------------------------------------------
# Copyright 2014 by Lilia Maliar, Serguei Maliar and Rafael Valero. All 
# rights reserved. The code may be used, modified and redistributed under  
# the terms provided in the file "License_Agreement.txt".
# -------------------------------------------------------------------------

function Smolyak_Elem_Anisotrop(Smol_elem_iso, smolyak_d, vector_mus_dimensions)


    points_dimensions = Int.(zeros(smolyak_d));   # This vector will tell how many
                                            # unidimensional elements in
                                            # each dimension we consider

    for i = 1:smolyak_d
        aux = vector_mus_dimensions[i];
        if aux == 0                            # If the approximation level in
                                                # the i-th dimension is 0, ...
            points_dimensions[i] = 1;         # The number of unidimensional 
                                                # elements is 1
        else                                    # If the approximation level in
                                                # the i-th dimension is not 0,...
            points_dimensions[i] = 2^(aux)+1; # Compute the number of unidimensional
                                                # elements using the formula
        end
    end

    for i = 1:smolyak_d
        aux1 = Smol_elem_iso[:,i] .<= points_dimensions[i];
        # If a subindex (i.e., number of elements) of the isotropic case is 
        # no larger than that of the anisotropic case, set an auxiliary variable
        # "aux1" to 1
        Smol_elem_iso = Smol_elem_iso[aux1, :];  
        # Keep only the vector of subindices for which aux=1
    end 
    return(Smol_elem_iso) 
end


# Smolyak_Grid is a routine that constructs a multidimensional Smolyak  
# grid in the hypercube [-1,1]^d; see "Smolyak method for solving dynamic 
# economic models: Lagrange interpolation, anisotropic grid and adaptive 
# domain" by Kenneth L. Judd, Lilia Maliar, Serguei Maliar and Rafael Valero,
# (2014), Journal of Economic Dynamics and Control 44, 92�123 (henceforth, 
# JMMV (2014)), Section 2.2.3 
#
# This version: Novenber 5, 2014. First version: December 17, 2012.
# -------------------------------------------------------------------------
# Input:   "d"         is the number of dimensions (the number of state 
#                      variables)
#          "mu"        is the level of approximation (in the anisotropic
#                      case, this is maximum level of approximation)
#          "smol_elem" is the matrix of the subindices of the Smolyak
#                      unidimensional elements; these elements can be either 
#                      isotropic (produced by Smolyak_Elem_Isotrop.m) or 
#                      anisotropic (produced by Smolyak_Elem_Anisotrop.m); 
#                      in the former case, the indices i1,...,id that jointly 
#                      satisfy the Smolyak rule, d<=|i|<=d+mu, where 
#                      |i|=i1+i2+...+id; see JMMV (2014), Section 3.2.3; 
#                      in the later case, they are a subset of the above 
#                      indices 
#
# Output:  "smol_grid" is the multidimensional Smolyak grid 
# -------------------------------------------------------------------------
# Copyright � 2014 by Lilia Maliar, Serguei Maliar and Rafael Valero. All 
# rights reserved. The code may be used, modified and redistributed under  
# the terms provided in the file "License_Agreement.txt".
# -------------------------------------------------------------------------


function Smolyak_Grid(d,mu,Smol_elem)

    # 1. Compute the vector of extrema of Chebyshev polynomials corresponding 
    # to the given level of Smolyak approximation mu
    # -----------------------------------------------------------------------

    # These points will be ordered as in Section 2.2.1 of JMMV(2014); e.g., for
    # mu=1, the set of points is {0,-1,1}
                
    points_1d = [];                    # Initially, the set of unidimensional 
                                        # points "points_1d" is empty; see JMMV  
                                        # (2014), Section 2.2.1
    i_max = mu+1;                      # The maximum subindex of unidimensional
                                        # set A_i whose points are used to
                                        # construct Smolyak grid of the given mu; 
                                        #  e.g., for mu=1, we consider up to 
                                        # A_i_max={-1,1} where i_max=1+1=2
    for i = 1:i_max                   # A subindex of a unidimensional set of
                                    # points                                
        # Compute the number of elements, m(i),(using m(i)=2^(i-1)+1) in the  
        # i-th unidimensional set of points; see Section 2.2.1 in JMMV (2014)
        #---------------------------------------------------------------------
        m_i = 1;                    # If i=1, then m(i)=1
        if i > 1
            m_i =  2^(i-1) + 1;     # If i>1, then m(i) = 2^(i-1)+1
        end                        
        
        # Construct the extrema of Chebyshev polynomials used as unidimensional 
        # grid points in the Smolyak method
        #---------------------------------------------------------------------
        if (m_i==1)
            extrem_Cheb_1d = 0;
        else 
        
            j = 1:m_i;                       
                                        # For j=1,...,m_i,...
            extrem_Cheb_1d = @. -cos(pi*(j-1)/(m_i-1));    
                                        # Chebyshev polynomials are defined in 
                                        # the interval [-1,1]
            extrem_Cheb_1d[abs.(extrem_Cheb_1d).<1e-12] .= 0.0;      
                                        # Round "extrem_Cheb_1d" to 0 if its  
                                        # absolute value is smaller than 1d-12
            extrem_Cheb_1d[abs.(1.0 .- extrem_Cheb_1d).<1e-12] .= 1.0;         
                                        # Round "extrem_Cheb_1d" to 1 if    
                                        # 1-extrem_Cheb_1d is smaller than 1d-12
            extrem_Cheb_1d[(1.0 .+extrem_Cheb_1d) .< 1e-12] .= -1.0;        
                                        # Round "extrem_Cheb_1d" to -1 if   
                                        # 1+extrem_Cheb_1d is smaller than 1d-12
                                        
        end 
    
        points_1d = vcat(points_1d, extrem_Cheb_1d);
                                    # Add to the previous set "points_1d" new 
                                    # points (given by extrema of unidimensional  
                                    # Chebyshev polynomials) as i increases
    
        points_1d = unique(points_1d,dims = 1);
                                    # Choose the unrepeated points and order 
                                    # them as in Section 2.2.1 of JMMV (2014); 
                                    
    # !!! NOTICE: For the versions of MATLAB older than MATLAB R2012a,  
    # option 'stable' might not work
    
    end              

        
    # 2. Construct the matrix multidimensional Smolyak grid points for the   
    # required level of Smolyak approximation, mu; see JMMV (2014), Sections 2.2.3 
    # for examples
    # -------------------------------------------------------------------------
    Smol_grid = zeros(size(Smol_elem));
                                # Initialize the matrix of multidimensional
                                # Smolyak grid points
                                
    numb_points_md = size(Smol_grid,1);
                                # Compute the number of multidimensional 
                                # Smolyak grid points                                 
    
    for jp = 1:numb_points_md     # For each multidimensional grid point, ...
                                
        index_row = Smol_elem[jp,:];
                                # Identify the subindex of the unidimensional
                                # grid point jp; this is a jp-th row of matrix 
                                # "Smol_elem"; 1-by-d
    
        for jd = 1:d              # For each dimension (state variable), ...
            n = index_row[jd];    # A subindex of a unidimensional grid point 
                                # in a dimension jd is denoted n
        
            Smol_grid[jp,jd] = points_1d[n];
                                # Find the corresponding unidimensional grid
                                # point in the vector "points_1d"
        end 
    end
    return(Smol_grid)
end


# Smolyak_Polynomial.m is a routine that constructs the multidimensional 
# basis functions of Smolyak polynomial of the approximation level that  
# corresponds to the previously constructed Smolyak (multidimensional) grid 
# points; see "Smolyak method for solving dynamic economic models: Lagrange 
# interpolation, anisotropic grid and adaptive domain" by Kenneth L. Judd, 
# Lilia Maliar, Serguei Maliar and Rafael, (2014). Journal of Economic 
# Dynamics and Control 44, 92�123 (henceforth, JMMV (2014)). 
#
# This version: Novenber 5, 2014. First version: May 30, 2011.
# -------------------------------------------------------------------------
# Inputs:  "points"    is the matrix of points in which the polynomial basis 
#                      functions must be evaluated; numb_pts-by-d
#          "d"         is the number of dimensions (state variables)
#          "smol_elem" is the matrix of the subindices of the Smolyak
#                      unidimensional elements; these elements can be either 
#                      isotropic (produced by Smolyak_Elem_Isotrop.m) or 
#                      anisotropic (produced by Smolyak_Elem_Anisotrop.m); 
#                      in the former case, the indices i1,...,id that jointly 
#                      satisfy the Smolyak rule, d<=|i|<=d+mu, where 
#                      |i|=i1+i2+...+id; see JMMV (2014), Section 3.2.3; 
#                      in the later case, they are a subset of the above 
#                      indices 
#
# Output:  "Smol_bases" is the matrix of multidimensional basis functions of 
#                      Smolyak polynomial of the given level of approximation, 
#                      evaluated in data matrix "points"
# -------------------------------------------------------------------------
# Copyright 2014 by Lilia Maliar, Serguei Maliar and Rafael Valero. All 
# rights reserved. The code may be used, modified and redistributed under  
# the terms provided in the file "License_Agreement.txt".
# -------------------------------------------------------------------------

function Smolyak_Polynomial(points,d,mu,Smol_elem)
 
    # Smolyak polynomial is given by the sum of multidimensional basis functions, 
    # multiplied by the coefficients; see formula (15) in JMMV (2014). By 
    # convention, the first basis function is given by 1 (unity). 

    # Unidimensional basis functions are given by Chebyshev polynomial bases; 
    # in JMMV (2014), a unidimensional Chebyshev polynomial basis function of   
    # degree n-1 is denoted by "phi_n", i.e., has a subindex n and we follow  
    # this notation here

    # 1. Construct the unidimensional basis functions and evaluate them in
    # all the points of matrix "points"
    # -------------------------------------------------------------------------
    i_max = mu+1;                   # The maximum subindex of unidimensional
                                    # set S_i whose points are used to
                                    # construct Smolyak grid; e.g., for mu=1, 
                                    # we consider up to S_i_max={0,-1,1} where
                                    # i_max=mu+1=1+1=2
                                    
    # Compute the number of elements in the i_max-th unidimensional set of 
    # elements S_i_max; this coincides with a maximum subindex of elements 
    # (unidimensional grid point or unidimensional basis function); 
    # see Section 2.2.1 in JMMV (2014)
    m_i_max = 1;
              # If i_max=1, then m(i_max)=1, i.e., set  
                                    # S_1={0} and the maximum subindex is 1
    if i_max > 1
        m_i_max =  2^(i_max-1) + 1;
                                    # If i_max>1, then m(i_max)= 2^(i_max-1)+1;
                                    # e.g., for S_2={0,-1,1}, the maximum
                                    # subindex is 3
    end
                                        
    numb_pts = size(points,1);      # Compute the number of points (rows),   
                                    # "numb_pts" in the matrix of points 
                                    # "points", in which the polynomial bases  
                                    # must be evaluated              
    phi = ones(numb_pts,d,m_i_max); 
                                    # Allocate memory to a matrix of the
                                    # unidimensional polynomial bases "phi_n",  
                                    # evaluated in all the points 
                            
    # For a polynomial bases "phi_n" with n=1, we have phi_n(x)=1 for all x; 
    # our phi(:,:,1) is a matrix of ones of size numb_pts-by-d by the above 
    # construction
                                    
    phi[:,:,2] .= points;            # For a polynomial bases "phi_n" with n=2, 
                                    # we have phi_n(x) is x; evaluating it in 
                                    # all the points gives us matrix "points"; 
                                    # numb_pts-by-d                    
    for j = 3:m_i_max              # For polynomial bases "phi_n", from n=3 to
                                # n=m_i_max, ...
        phi[:,:,j] = 2 .* phi[:,:,2].*phi[:,:,j-1] .- phi[:,:,j-2]; 
                                # Use the recurrence formula to compute the 
                                # Chebyshev polynomial basis functions of 
                                # the degrees from 2 to m_i_max-1
    end 

        
    # 2. Form the multidimensional polynomial bases of Smolyak polynomial of the 
    # required level of Smolyak approximation; see JMMV (2014), Sections 3.3.3 
    # and 3.4.2 for examples
    # ----------------------------------------------------------------------
    Smol_bases = [];              # Initially, the matrix of multidimensional 
                                # polynomial bases is empty
                                
    numb_terms = size(Smol_elem,1);
                                # Compute the number of terms (i.e., multi-
                                # dimensional polynomial bases) in Smolyak 
                                # polynomial                                
    
    for jt = 1:numb_terms         # For each term of Smolyak polynomial, ...
                                
        index_row = Smol_elem[jt,:];
                                # Identify the subindices of the unidimensional
                                # basis function that constitute a jt-th multi-
                                # dimensional basis function; this is a jt-th
                                # row of matrix "Smol_elem"; 1-by-d
        product = ones(numb_pts,1); 
                                # Initialize a vector of products of unidimen-
                                # sional basis functions; numb_pts-by-1
        for jd = 1:d              # For each dimension (state variable), ...
            n = index_row[jd];    # A subindex of a unidimensional basis function 
                                # phi_n in a dimension jd is denoted n
            if n ≠ 1;            # If the subindex of unidimensional basis 
                                # function is not equal to unity, ...
                product = product.*phi[:,jd,n];
                                # Compute the product of basis functions
                                
                                # Otherwise, i.e., if n = 1, there is no need to compute the 
                                # product of unidimensional basis functions, as it's equal to unity 
            end
        end
        
        if isempty(Smol_bases)
            Smol_bases = product' 
        else
            Smol_bases = vcat(Smol_bases,product');
                                    # Attach to the previously obtained matrix of 
                                    # multidimensional basis functions a new
                                    # product of unidimensional basis functions;
                                    # e.g., for mu=1 and d=2, basis_bs is of
                                    # size numb_pts-by-5
        end
    end
    return(Smol_bases')
end

function  Smolyak_Polynomial2(points,d,numb_pts,numb_terms,mu,smol_elem)
    smol_bases = ones(numb_pts,numb_terms)  
    # numb_pts = size(points,1)       # Compute the number of points (rows),   
                                    # "numb_pts" in the matrix of points 
                                    # "points", in which the polynomial bases  
                                    # must be evaluated              
    # allocate(phi(1:numb_pts,1:d,1:m_i_max))
    # Smolyak polynomial is given by the sum of multidimensional basis functions, 
    # multiplied by the coefficients; see formula (15) in JMMV (2014). By 
    # convention, the first basis function is given by 1 (unity). 

    # Unidimensional basis functions are given by Chebyshev polynomial bases; 
    # in JMMV (2014), a unidimensional Chebyshev polynomial basis function of   
    # degree n-1 is denoted by "phi_n", i.e., has a subindex n and we follow  
    # this notation here

    # 1. Construct the unidimensional basis functions and evaluate them in
    # all the points of matrix "points"
    # -------------------------------------------------------------------------
    i_max = mu+1                    # The maximum subindex of unidimensional
                                    # set S_i whose points are used to
                                    # construct Smolyak grid; e.g., for mu=1, 
                                    # we consider up to S_i_max={0,-1,1} where
                                    # i_max=mu+1=1+1=2
                                    
    # # Compute the number of elements in the i_max-th unidimensional set of 
    # # elements S_i_max; this coincides with a maximum subindex of elements 
    # # (unidimensional grid point or unidimensional basis function); 
    # # see Section 2.2.1 in JMMV (2014)

    m_i_max =  2^(i_max-1) + 1
    if (i_max == 1)  
        m_i_max = 1         # If i_max=1,  m(i_max)=1, i.e., set  
                            # S_1={0} and the maximum subindex is 1
    end                   # If i_max>1,  m(i_max)= 2^(i_max-1)+1;
                            # e.g., for S_2={0,-1,1}, the maximum
                            # subindex is 3

                                    
    phi = ones(numb_pts, d, 2^mu+1);
                                    # Allocate memory to a matrix of the
                                    # unidimensional polynomial bases "phi_n",  
                                    # evaluated in all the points 
                        
    # For a polynomial bases "phi_n" with n=1, we have phi_n(x)=1 for all x; 
    # our phi(:,:,1) is a matrix of ones of size numb_pts-by-d by the above 
    # construction
                                
    phi[:,:,2] .= points;            # For a polynomial bases "phi_n" with n=2, 
                                    # we have phi_n(x) is x; evaluating it in 
                                    # all the points gives us matrix "points"; 
                                    # numb_pts-by-d                    
    for j = 3:m_i_max                # For polynomial bases "phi_n", from n=3 to
                                    # n=m_i_max, ...
        phi[:,:,j] = 2 .* phi[:,:,2] .* phi[:,:,j-1] .- phi[:,:,j-2];
                                    # Use the recurrence formula to compute the 
                                    # Chebyshev polynomial basis functions of 
                                    # the degrees from 2 to m_i_max-1
    end 

    
    # 2. Form the multidimensional polynomial bases of Smolyak polynomial of the 
    # required level of Smolyak approximation; see JMMV (2014), Sections 3.3.3 
    # and 3.4.2 for examples
    # ----------------------------------------------------------------------
    # Smol_bases = [];             # Initially, the matrix of multidimensional 
                                # polynomial bases is empty
                                
    # numb_terms = size(smol_elem,1);
                                # Compute the number of terms (i.e., multi-
                                # dimensional polynomial bases) in Smolyak 
                                # polynomial                                

    # allocate(pproduct(1:numb_pts))
    
    for jt = 1:numb_terms          # For each term of Smolyak polynomial, ...
                                    
        index_row = smol_elem[jt,:]
                                # Identify the subindices of the unidimensional
                                # basis function that constitute a jt-th multi-
                                # dimensional basis function; this is a jt-th
                                # row of matrix "smol_elem"; 1-by-d
        pproduct = ones(numb_pts);
                                # Initialize a vector of pproducts of unidimen-
                                # sional basis functions; numb_pts-by-1
        for jd = 1:d              # For each dimension (state variable), ...
            n = index_row[jd]    # A subindex of a unidimensional basis function 
                                # phi_n in a dimension jd is denoted n
            if n ≠ 1              # If the subindex of unidimensional basis 
                                # function is not equal to unity, ...
                pproduct = pproduct .* phi[:,jd,n];
                                # Compute the pproduct of basis functions
                                
            # Otherwise, i.e., if n = 1, there is no need to compute the 
            # pproduct of unidimensional basis functions, as it's equal to unity 
            end
        end

        smol_bases[:,jt] = pproduct

    #  if (allocated(Smol_bases)) 
    #  Smol_bases_temp = Smol_bases
    #  deallocate(smol_bases)
    #  allocate(Smol_bases(1:numb_pts,size(Smol_bases_temp,2)+1))
    #  Smol_bases(:,1:size(Smol_bases_temp,2)) = Smol_bases_temp
    #  Smol_bases(:,size(Smol_bases_temp,2)+1) = pproduct
    #  else 
    #  allocate(Smol_bases(1:numb_pts,1:1))
    #  Smol_bases(1:numb_pts,1) = pproduct
    #  end
                                # Attach to the previously obtained matrix of 
                                # multidimensional basis functions a new
                                # product of unidimensional basis functions;
                                # e.g., for mu=1 and d=2, basis_bs is of
                                # size numb_pts-by-5
    end
    return(smol_bases)
end 

