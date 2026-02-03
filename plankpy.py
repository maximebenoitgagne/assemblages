#!/usr/bin/env python3

"""Classes to describe the output of a simulation by MITgcm

The behaviour of the classes is defined if the configuration of MITgcm
was a 1-D configuration similar to the configuration used for
https://github.com/maximebenoitgagne/timing (Benoît-Gagné et al., 2024),
DOI: https://doi.org/10.5281/zenodo.10642024.
"""

# author: Maxime Benoit-Gagne - ULaval - Canada
# date of creation: February 21, 2024
#
# Python from Anaconda.
#
# % python
# Python 3.8.18 | packaged by conda-forge | (default, Dec 23 2023, 17:25:47) 
# [Clang 16.0.6 ] on darwin
# Type "help", "copyright", "credits" or "license" for more information.

########### importing modules ###########
import numpy as np

import netcdf_tools
import vstats_pos

########### constants ###########

__FIRST_TRACER_PLANKTON__=21
__MAP_CXX_GROUP_SIZE__={
    1:('pro','0.6'),
    2:('syn','0.9'),
    3:('pico','1.4'),
    4:('pico','2.0'),
    5:('other','3.0'),
    6:('other','4.5'),
    7:('other','6.6'),
    8:('other','10'),
    9:('other','15'),
    10:('diazo','3.0'),
    11:('diazo','4.5'),
    12:('diazo','6.6'),
    13:('diazo','10'),
    14:('tricho','15'),
    15:('diatoms','6.6'),
    16:('diatoms','10'),
    17:('diatoms','15'),
    18:('diatoms','22'),
    19:('diatoms','32'),
    20:('diatoms','47'),
    21:('diatoms','70'),
    22:('diatoms','104'),
    23:('diatoms','154'),
    24:('dino','6.6'),
    25:('dino','10'),
    26:('dino','15'),
    27:('dino','22'),
    28:('dino','32'),
    29:('dino','47'),
    30:('dino','70'),
    31:('dino','104'),
    32:('dino','154'),
    33:('dino','228'),
    34:('zoo','6.6'),
    35:('zoo','10'),
    36:('zoo','15'),
    37:('zoo','22'),
    38:('zoo','32'),
    39:('zoo','47'),
    40:('zoo','70'),
    41:('zoo','104'),
    42:('zoo','154'),
    43:('zoo','228'),
    44:('zoo','338'),
    45:('zoo','502'),
    46:('zoo','744'),
    47:('zoo','1103'),
    48:('zoo','1636'),
    49:('zoo','2425'),  
}
__MOLAR_MASS_C__=12.0107 # g C (mol C)^-1

########### class ###########
class Ecosystem(object):
    """
    Information about an ecosystem in the simulation.

    Attributes:
        array1d_ip_p(array-like):
            Array of 1 dimension.
            The first dimension is the indices of the plankton type
            from 0 to 48.
            They correspond to cxx-1.
            The values are each plankton type.
    """
    def __init__(self,
                 carfile,
                 gridfile,
                 array1d_iT1y_iT,
                 depth_end):
        """
        The __init__ method.

        Args:
            carfile(str):
                The netCDF4 file for biomass from a simulation by the
                biogeochemical component of MITgcm.
            gridfile(str):
                The netCDF4 file for the grid of a simulation by
                MITgcm.
                The file shall include a variable RC for the
                coordinates of the cell center.
            array1d_iT1y_iT(array-like):
                Array of 1 dimension.
                The first dimension (of length 365) is the indices of
                the indices of the days.
                The values are the indices of the days.
                For example, if carfile contains daily data over 10
                years and we want the tenth year, we need to get
                the values at the indices 3285 to 3649.
                array1d_iT1y_iT is that 3285 to 3649.
            depth_end(float):
                Depth to which integration is made.
        """
        nplankton=len(__MAP_CXX_GROUP_SIZE__)
        self.array1d_ip_p=np.zeros(
            shape=nplankton,
            dtype=object)
        for ip in range(0,nplankton):
            cxx=ip+1
            group_name,size_name=Plankton.get_group_size(cxx)
            plankton=Plankton(carfile,
                              gridfile,
                              array1d_iT1y_iT,
                              depth_end,
                              group_name,
                              size_name)
            self.array1d_ip_p[ip]=plankton

    def get_array2d_idepth_iT_c(self,array1d_ip_gs):
        """
        Return the sum of the biomass concentrations of each plankton.
        
        Args:
            array1d_ip_gs(array-like):
                Array of 1 dimension.
                The first dimension is the indices of plankton
                The values are tuples containing the group and the size
                class of the plankton.
                
        Returns:
            array-like(int):
            Array of 2 dimensions.
            The first dimension is the indices of the depths.
            The second dimension is the indices of the time depths.
            The values are the sum of the biomass concentrations of each
            plankton (mg C m^-3).
                
        Raises:
            ValueError: If a group_name and a size_name do not
            correspond to a valid plankton type.
        """
        assert len(self.array1d_ip_p)>0, \
        'The Ecosystem contains no plankton.'
        p=self.array1d_ip_p[0]
        ndepths,ndays=p.array2d_idepth_iT_c.shape
        array2d_idepth_iT_c=np.zeros((ndepths,ndays))
        for ip,group_size in enumerate(array1d_ip_gs):
            p=self.get_plankton(*group_size)
            array2d_idepth_iT_ctempo=p.array2d_idepth_iT_c
            array2d_idepth_iT_c= \
            array2d_idepth_iT_c+array2d_idepth_iT_ctempo
        return array2d_idepth_iT_c

    def get_group(self,group_name):
        """
        Return the list of the plankton in the group.

        Args:
            group_name(str):
                Short name of the plankton group.
                Accepted values are 'pro', 'syn', 'pico', 'other',
                'diazo', 'tricho', 'diatoms', 'dino' and 'zoo'.

        Returns:
            list(Plankton):
            The list of the plankton in the group.
            
        Raises:
            ValueError: If the group_name does not correspond to a valid 
            group of plankton.
        """
        keys_group_name=[k for k, v in __MAP_CXX_GROUP_SIZE__.items()
                         if v[0]==group_name]
        n=len(keys_group_name)
        if n==0:
            errmsg="""The accepted values of group_name are:"""
            prev_str=""
            for k,v in sorted(__MAP_CXX_GROUP_SIZE__.items()):
                str="\ngroup_name="+"{0:<9}".format("'"+v[0]+"',")
                if str!=prev_str:
                    errmsg=errmsg+str
                prev_str=str
            errmsg=errmsg+"\nGot group_name={0}".format(group_name)
            raise ValueError(errmsg)
        group=[]
        for cxx in keys_group_name:
            group_name,size_name=Plankton.get_group_size(cxx)
            plankton=self.get_plankton(group_name,size_name)
            group.append(plankton)
            assert isinstance(plankton,Plankton), \
            """Plankton group_name={0} and size_name={1} exists in
the dictionary of module plankpy but was not initialized in the
ecosystem.""".format(group_name,size_name)
        return group

    def get_group_sumc(self,group_name):
        """
        Return the vertically integrated group biomass.

        Args:
            group_name(str):
                Short name of the plankton group.
                Accepted values are 'pro', 'syn', 'pico', 'other',
                'diazo', 'tricho', 'diatoms', 'dino' and 'zoo'.

        Returns:
            (array-like):
            Array of 1 dimension.
            The first dimension is the indices of the time steps.
            The values are the vertically integrated group biomass
            (0-100m, in mg C m^-2).
            
        Raises:
            ValueError: If the group_name does not correspond to a valid 
            group of plankton.
        """
        group=self.get_group(group_name)
        assert len(group)>1,"""The accepted values of group_name are
        'pro', 'syn', 'pico', 'other','diazo', 'tricho', 'diatoms', 'dino'
        and 'zoo'. Got group_name={0}.format(group_name)"""
        ndays=len(group[0].array1d_iT_sumc)
        group_sumc=np.zeros(ndays)
        for i in range(len(group)):
            group_sumc=group_sumc+group[i].array1d_iT_sumc
        return group_sumc
    
    def get_plankton(self,group_name,size_name):
        """
        Return the corresponding plankton.
        
        Args:
            group_name(str):
                Short name of the plankton group.
                Accepted values are 'pro', 'syn', 'pico', 'other',
                'diazo', 'tricho', 'diatoms', 'dino' and 'zoo'.
            size_name(str):
                Class size of the plankton.
                Accepted values are '0.6', '0.9', '1.4', '2.0', '3.0',
                '4.5', '6.6', '10', '15', '22', '32', '47', '70',
                '104', '154', '228', '338', '502', '744', '1103',
                '1636', '2425'.
                
        Returns:
            Plankton: The corresponding plankton.
            
        Raises:
            ValueError: If the group_name and the size_name do not
            correspond to a valid plankton type.
        """
        cxx=Plankton.get_cxx(group_name,size_name) # 1-based
        ip=cxx-1 # 0-based
        plankton=self.array1d_ip_p[ip]
        assert isinstance(plankton,Plankton), \
        """Plankton group_name={0} and size_name={1} exists in
the dictionary of module plankpy but was not initialized in the
ecosystem."""\
        .format(group_name,size_name)
        return plankton

    def __repr__(self):
        """
        Compute the "official" string representation of an object.

        Returns:
            str: The "official" string representation of an object.
        """
        return """Ecosystem(array1d_ip_p={0.array1d_ip_p!r})""".format(self)
    
    def __str__(self):
        """
        Compute the "informal" string representation of an object.

        Returns:
            str:
            The "informal" string representation of an object.
        """
        s=repr(self)
        s=s.replace(self.__class__.__name__,"")
        return s

########### class ###########
class Plankton(object):
    """
    Information about one plankton type in the simulation.

    Attributes:
        array1d_iT_sumc (array-like):
            Array of 1 dimension.
            The first dimension is the indices of the time steps.
            The values are the vertically integrated plankton biomass
            (0-100m, in mg C m^-2).
        array2d_idepth_iT_c (array-like)
            Array of 2 dimensions.
            The first dimension is the indices of the depths.
            The second dimension is the indices of the time depths.
            The values are the plankton biomass concentration
            (mg C m^-3).
    """

    def __init__(self,
                 carfile,
                 gridfile,
                 array1d_iT1y_iT,
                 depth_end,
                 group_name,
                 size_name):
        """
        The __init__ method.

        Args:
            carfile(str):
                The netCDF4 file for biomass from a simulation by the
                biogeochemical component of MITgcm.
            gridfile(str):
                The netCDF4 file for the grid of a simulation by
                MITgcm.
                The file shall include a variable RC for the
                coordinates of the cell center.
            array1d_iT1y_iT(array-like):
                Array of 1 dimension.
                The first dimension (of length 365) is the indices of
                the indices of the days.
                The values are the indices of the days.
                For example, if carfile contains daily data over 10
                years and we want the tenth year, we need to get
                the values at the indices 3285 to 3649.
                array1d_iT1y_iT is that 3285 to 3649.
            depth_end(float):
                Depth to which integration is made.
            group_name(str):
                Short name of the plankton group.
                Accepted values are 'pro', 'syn', 'pico', 'other',
                'diazo', 'tricho', 'diatoms', 'dino' and 'zoo'.
            size_name(str):
                Class size of the plankton.
                Accepted values are '0.6', '0.9', '1.4', '2.0', '3.0',
                '4.5', '6.6', '10', '15', '22', '32', '47', '70',
                '104', '154', '228', '338', '502', '744', '1103',
                '1636', '2425'.
        """
        tracer_name=Plankton.get_tracer_name(group_name,size_name)
        array2d_idepth_iT_c\
        =netcdf_tools.read_netcdf(carfile,tracer_name)\
        .squeeze().transpose()
        array2d_idepth_iT_c=array2d_idepth_iT_c[:,array1d_iT1y_iT]
        array2d_idepth_iT_c[-1,:]=np.nan
        self.array2d_idepth_iT_c=array2d_idepth_iT_c*__MOLAR_MASS_C__
        drF=netcdf_tools.read_netcdf(gridfile,'drF')
        self.array1d_iT_sumc=vstats_pos.vint(
            array2d_idepth_iT_tracer=self.array2d_idepth_iT_c,
            array1d_idepth_delR=drF,
            depth_end=depth_end)
    
    @staticmethod
    def get_cxx(group_name,size_name):
        """
        Return the number of the plankton.
        
        Return the last two digits in the name of the variable
        corresponding to the group and size of the plankton.
        For example, if the plankton corresponds to variable c15,
        the return value will be 15.
        
        Args:
            group_name(str):
                Short name of the plankton group.
                Accepted values are 'pro', 'syn', 'pico', 'other',
                'diazo', 'tricho', 'diatoms', 'dino' and 'zoo'.
            size_name(str):
                Class size of the plankton.
                Accepted values are '0.6', '0.9', '1.4', '2.0', '3.0',
                '4.5', '6.6', '10', '15', '22', '32', '47', '70',
                '104', '154', '228', '338', '502', '744', '1103',
                '1636', '2425'.
                
        Returns:
            int: The number of the plankton.
            -21 if the group_name and the size_name do not
            correspond to a valid plankton type.
            
        Raises:
            ValueError: If the group_name and the size_name do not
            correspond to a valid plankton type.
        """
        cxx=-21
        group_size=(group_name,size_name)
        keys=list(__MAP_CXX_GROUP_SIZE__.keys())
        values=list(__MAP_CXX_GROUP_SIZE__.values())
        if group_size in values:
            cxx=keys[values.index(group_size)]
        else:
            errmsg="""The accepted values of \
group_name, size_name are:"""
            for k,v in sorted(__MAP_CXX_GROUP_SIZE__.items()):
                str="\ngroup_name="+"{0:<9}".format("'"+v[0]+"'")\
                +" and size_name="+"{0:<9}".format("'"+v[1]+"',")
                errmsg=errmsg+str
            errmsg=errmsg+"\nGot group_name={0} and size_name={1}"\
                    .format(group_name,size_name)
            raise ValueError(errmsg)
        return cxx
    
    @staticmethod
    def get_group_size(cxx):
        """
        Return the group and the size of the plankton number cxx.
        
        Args:
            cxx(int):
                The number of the plankton.
                
        Returns:
            (str,str): The group and size class of plankton cxx.

        Raises:
            ValueError: If cxx is not a valid number of plankton.
        """
        if cxx not in __MAP_CXX_GROUP_SIZE__:
            errmsg="cxx={0} is not a valid plankton number".format(cxx)
            raise ValueError(errmsg)
        return __MAP_CXX_GROUP_SIZE__[cxx]

    @staticmethod
    def get_isize_per_group(cxx):
        """
        Return the index of the size class relative to the group
        of the plankton number cxx.
                
        Args:
            cxx(int):
                The number of the plankton.
                    
        Returns:
            int: returns the index (0-based) of the size class
            relative to the group of the plankton number cxx.
            -1 if cxx is not a valid number of plankton.
                    
        Raises:
            ValueError: If cxx is not a valid number of plankton.
        """
        group_name,size_name=Plankton.get_group_size(cxx)
        keys_group_name=[k for k, v in __MAP_CXX_GROUP_SIZE__.items()
                         if v[0]==group_name]
        cxx_first=keys_group_name[0]
        return cxx-cxx_first

    @staticmethod
    def get_nsize_per_group(group_name):
        """
        Return the number of size classes for the group.
            
        Args:
            group_name(str):
                Short name of the plankton group.
                Accepted values are 'pro', 'syn', 'pico', 'other',
                'diazo', 'tricho', 'diatoms', 'dino' and 'zoo'.
                
        Returns:
            int: Return the number of size classes for the group.
            -1 if the group_name does not correspond to a valid group.
                
            Raises:
                ValueError: If the group_name does not correspond
                to a valid group.
        """
        keys_group_name=[k for k, v in __MAP_CXX_GROUP_SIZE__.items()
                         if v[0]==group_name]
        n=len(keys_group_name)
        if n==0:
            n=-1
            errmsg="""The accepted values of group_name are:"""
            prev_str=""
            for k,v in sorted(__MAP_CXX_GROUP_SIZE__.items()):
                str="\ngroup_name="+"{0:<9}".format("'"+v[0]+"',")
                if str!=prev_str:
                    errmsg=errmsg+str
                prev_str=str
            errmsg=errmsg+"\nGot group_name={0}".format(group_name)
            raise ValueError(errmsg)
        return n
    
    @staticmethod
    def get_tracer_name(group_name,size_name):
        """
        Return the tracer name corresponding to the group and size
        of the plankton.
        
        Args:
            group_name(str):
                Short name of the plankton group.
                Accepted values are 'pro', 'syn', 'pico', 'other',
                'diazo', 'tricho', 'diatoms', 'dino' and 'zoo'.
            size_name(str):
                Class size of the plankton.
                Accepted values are '0.6', '0.9', '1.4', '2.0', '3.0',
                '4.5', '6.6', '10', '15', '22', '32', '47', '70',
                '104', '154', '228', '338', '502', '744', '1103',
                '1636', '2425'.
                
        Returns:
            str: The tracer name TRACxx.
        """
        cxx=Plankton.get_cxx(group_name,size_name)
        tracer_number=cxx+__FIRST_TRACER_PLANKTON__-1
        tracer_name="TRAC{0:02}".format(tracer_number)
        return tracer_name

    def __repr__(self):
        """
        Compute the "official" string representation of an object.

        Returns:
            str: The "official" string representation of an object.
        """
        return """Plankton(array2d_idepth_iT_c={0.array2d_idepth_iT_c!r},
array1d_iT_sumc={0.array1d_iT_sumc!r})""".format(self)
    
    def __str__(self):
        """
        Compute the "informal" string representation of an object.

        Returns:
            str:
            The "informal" string representation of an object.
        """
        s=repr(self)
        s=s.replace(self.__class__.__name__,"")
        return s
