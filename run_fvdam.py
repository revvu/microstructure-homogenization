# Reevu Adakroy
# July 28, 2022

import shutil
import matlab.engine
import os

LOP_naming = ['SIGMA11','SIGMA22','SIGMA33','SIGMA23','SIGMA13','SIGMA12']

increment_count = 25

def input_filename(time_saved):
    return 'static/Generated Microstructures/FVDAM Input/generated_microstructure_input_'+time_saved+'.fgm'

def output_filename(time_saved, LOP):
    return 'static/Generated Microstructures/FVDAM Output/generated_microstructure_output_'+time_saved+'_LOP_'+str(LOP)+'.out'

def main(time_saved, LOP):
    os.chdir('./static/FVDAM')
    print('Engine starting...')
    eng = matlab.engine.start_matlab()
    print('Engine started.')
    os.chdir('..')
    os.chdir('..')

    try:
        shutil.copy(input_filename(time_saved), './static/FVDAM/generated_ruc_46_bal_insitu.fgm')
        print('file copied successfully')
        if LOP==1:   eng.fvdam_global_exec1(nargout=0)
        elif LOP==2: eng.fvdam_global_exec2(nargout=0)
        elif LOP==3: eng.fvdam_global_exec3(nargout=0)
        elif LOP==4: eng.fvdam_global_exec4(nargout=0)
        elif LOP==5: eng.fvdam_global_exec5(nargout=0)
        elif LOP==6: eng.fvdam_global_exec6(nargout=0)

        shutil.copy('./static/FVDAM/generated_ruc_46_bal_insitu.out', output_filename(time_saved,LOP))
    except FileNotFoundError: pass

    return output_filename(time_saved,LOP)

if __name__ == '__main__': main('1659071042',1)