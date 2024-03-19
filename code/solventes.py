from collections import defaultdict
import sys
import time
import os
import pandas as pd


class solventes(object):

    def get_no_of_atoms(self, **kwargs):

        path_file = kwargs.get("path_file", 3)

        f = open(path_file, "r")
        for line in f:
            if "No. of atoms" in line:
                N = line.split()[-1]
                break
        f.close()

        return int(N)
    

    def get_xyzs_from_nwchem(self, **kwargs):
        """
        return: pandas dataframe

        inputs:

        path_file
        atom_number
        
        """

        path_file = kwargs.get("path_file", "salida.out")
        atom_number = self.get_no_of_atoms(path_file=path_file)
    
        au2aFactor = 0.5291771057875306

        dic = defaultdict(dict)

        ij = 0
        f = open(path_file, "r")
        idx = 0
        for line in f:
            if "DFT ENERGY GRADIENTS" in line:
                line = next(f)
                line = next(f)
                line = next(f)
                for i in range(atom_number):
                    line = next(f)
                    sp = line.split()
                    
                    ida = sp[0]
                    atom = sp[1]
                    x = sp[2]
                    y = sp[3]
                    z = sp[4]
                    fx = sp[5]
                    fy = sp[6]
                    fz = sp[7]
                    frame = ij
                    
                    dic["ida"][idx] = ida
                    dic["atom"][idx] = atom
                    dic["x"][idx] = float(x)*au2aFactor
                    dic["y"][idx] = float(y)*au2aFactor
                    dic["z"][idx] = float(z)*au2aFactor
                    dic["fx"][idx] = float(fx)
                    dic["fy"][idx] = float(fy)
                    dic["fz"][idx] = float(fz)        
                    dic["frame"][idx] = int(frame)
                    idx+=1
                
                ij+=1
        
        f.close()

        return pd.DataFrame(dic)




    def get_xyz_movie(self, **kwargs):

        path_file = kwargs.get("path_file", "salida.out")
        path_movie = kwargs.get("path_movie", "movie.xyz")

        df = self.get_xyzs_from_nwchem(path_file=path_file)

        f = open(path_movie, "w")
        
        frames = df["frame"].unique()
        N = len(df[df["frame"] == 0])
        for frame in frames:
            
            tmp_df = df[df["frame"] == frame]
            f.write(f"{len(tmp_df)}")
            f.write(f"\n  frame   {frame}\n")
            for idx in tmp_df.index:
                
                atom = tmp_df["atom"].loc[idx]
                x = tmp_df["x"].loc[idx]
                y = tmp_df["y"].loc[idx]
                z = tmp_df["z"].loc[idx]
                
                f.write(f"{atom}")
                f.write(f"{x:10.4f}{y:10.4f}{z:10.4f}\n")
                
            
            
        f.close()


        return


    def get_energies_df(self, **kwargs):

        path_file = kwargs.get("path_file", "salida.out")
        label = kwargs.get("label", "solvente")
        
        txt_2find = "Total DFT energy"
        
        N = self.get_no_of_atoms(path_file=path_file)
        
        au2aFactor = 0.5291771057875306
        dic2 = defaultdict(dict)
        
        dfxyz_ = pd.DataFrame({})
        idy = 0 # frames
        f = open(path_file, 'r')
        
        for line in f:
            
            if txt_2find in line:
                energy = line.split()[4]
                dic2["energy"][idy] = float(energy)
            
            if "DFT ENERGY GRADIENTS" in line:
                line = next(f)
                line = next(f)
                line = next(f)
                
                dic = defaultdict(dict)
                
                idx=0
                for i in range(N):
                    line = next(f)
                    sp = line.split()
                    
                    ida = sp[0]
                    atom = sp[1]
                    x = sp[2]
                    y = sp[3]
                    z = sp[4]
                    fx = sp[5]
                    fy = sp[6]
                    fz = sp[7]
                    
                    dic["ida"][idx] = ida
                    dic["atom"][idx] = atom
                    dic["x"][idx] = float(x)*au2aFactor
                    dic["y"][idx] = float(y)*au2aFactor
                    dic["z"][idx] = float(z)*au2aFactor
                    dic["fx"][idx] = float(fx)
                    dic["fy"][idx] = float(fy)
                    dic["fz"][idx] = float(fz)        
                    dic["frame"][idx] = idy  
                    
                    idx+=1
                    
                dfxyz = pd.DataFrame(dic)
                dfxyz_ = pd.concat([dfxyz_, dfxyz], ignore_index=True)        
         
                idy+=1
                
                
        f.close()

        dfe = pd.DataFrame(dic2)


        return dfe,dfxyz_
        




