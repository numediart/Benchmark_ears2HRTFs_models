import torch
"""
Convert HRTF Models to SOFA Format.
This script loads pre-trained PyTorch models designed for HRTF estimation, 
generates Head-Related Transfer Functions (HRTFs) for a dataset of subjects, 
applies minimum phase reconstruction via Hilbert transform, adds pure delay 
based on geometric calculation, converts the frequency domain data to the 
time domain (HRIR), and saves the results as SOFA files.
Key Features:
- Iterates through a directory of trained model weights.
- Handles different model architectures (e.g., specific logic for 'Manlin-Zhao').
- Calculates pure geometric delay based on incidence angles.
- Reconstructs phase information using the Hilbert transform.
- Exports data using the `sofar` library.
Dependencies:
    torch, os, sofar, scipy, numpy, tqdm
    """

"""
    Helper function to extract epoch numbers from filenames for sorting.
    Assumes filenames contain epoch information around index 11-13.
    Args:
        name (str): The filename string.
    Returns:
        int: The extracted epoch number.
    """

# def pure_dealy_com(Incidence_angle, speed_sound=343, distance_head_center_source=1.5, distance_ear_head=0.09):
"""
    Computes the complex exponential component representing pure time delay due to geometry.
    Calculates the distance from source to ears based on a spherical head model approximation
    and converts these distances into phase shifts in the frequency domain.
    Args:
        Incidence_angle (np.ndarray): Array of shape (N, 2) containing [azimuth, elevation] in degrees.
        speed_sound (float, optional): Speed of sound in m/s. Defaults to 343.
        distance_head_center_source (float, optional): Distance from head center to source in meters. Defaults to 1.5.
        distance_ear_head (float, optional): Distance from ear to head center (radius) in meters. Defaults to 0.09.
    Returns:
        np.ndarray: A complex array of shape (9, N_angles, 2, 256) representing the phase shift factors 
                    associated with the calculated pure delays.
                    Note: The first dimension (9) appears to correspond to a batch or fixed repetition size 
                    specific to the calling context.
    """


# def HRTF_TO_HRIR_with_phase(HRTF_log_magnitude: np.ndarray, incidence_angle: np.ndarray, model_name: str):
"""
    Converts HRTF log-magnitude predictions to time-domain HRIRs and saves them as SOFA files.
    This process involves:
    1. Converting log-magnitude to linear amplitude.
    2. Mirroring the spectrum to create a symmetric Real-FFT format.
    3. Reconstructing phase using the Hilbert transform (minimum phase assumption).
    4. Applying pure geometric delay calculated via `pure_dealy_com`.
    5. Performing Inverse FFT to get the Impulse Response (HRIR).
    6. Packaging the data into a 'SimpleFreeFieldHRIR' SOFA object and saving it to disk.
    Args:
        HRTF_log_magnitude (np.ndarray): Tensor of shape (N_batch, 793, 2, 128) containing predicted log magnitudes.
        incidence_angle (np.ndarray): Array of angles corresponding to the HRTF positions.
        model_name (str): Name identifier for the model, used in the output filename.
    Returns:
        str: Returns a string literal "uwu" upon completion.
    """
from Datasets import HRTF_mesh_dataset
import os 
import sofar as sf
from scipy.signal import hilbert
import numpy as np
from tqdm import tqdm

def func_sort(name):
    str_epoch=name[11:13]
    if str_epoch[1]=='_':
        str_epoch=str_epoch[0]
    int_epoch= int(str_epoch)
    return int_epoch

def pure_dealy_com(Incidence_angle, speed_sound=343,distance_head_center_source=1.5,distance_ear_head=0.09):
    pure_delays=np.zeros((Incidence_angle.shape[0],2))
    pure_delays_exp_fin_form=np.zeros((9,Incidence_angle.shape[0],2,256),dtype=complex)
    omega=np.linspace(0,44100*2*np.pi,256)
    omega=np.matlib.repmat(omega,9,1)
    distance_delay=np.zeros((Incidence_angle.shape[0],2)) # azimuth,elevation (0,1)
    for i in range(Incidence_angle.shape[0]):
        distance_delay[i,0]=np.sqrt(distance_head_center_source**2 + distance_ear_head**2 - 2*distance_head_center_source*distance_ear_head*np.cos(np.pi-np.deg2rad(Incidence_angle[i,1]))*np.sin(np.deg2rad(Incidence_angle[i,0])))
        distance_delay[i,1]=np.sqrt(distance_head_center_source**2 + distance_ear_head**2 - 2*distance_head_center_source*distance_ear_head*np.cos(np.deg2rad(Incidence_angle[i,1]))*np.sin(np.deg2rad(Incidence_angle[i,0])))
        pure_delays[i,:]=distance_delay[i,:]/speed_sound
        pure_delays_exp_fin_form[:,i,0,:]=np.exp(1j*omega* pure_delays[i,0])
        pure_delays_exp_fin_form[:,i,1,:]=np.exp(1j*omega* pure_delays[i,1])
    return pure_delays_exp_fin_form 

def HRTF_TO_HRIR_with_phase(HRTF_log_magnitude:np.ndarray,incidence_angle:np.ndarray,model_name:str): #(N_batch,793,2,128)
    dist=np.ones((793,1))*1.5
    HRTF_fft_amplitude = 10 ** (HRTF_log_magnitude / 20)
    HRTF_fft_fin = np.concatenate((HRTF_fft_amplitude, np.flip(HRTF_fft_amplitude, axis=-1)), axis=-1)
    Phase_phi=-hilbert(np.log10(HRTF_fft_fin), N=None, axis=-1)
    HRTF_fft_fin = HRTF_fft_fin * np.exp(1j * Phase_phi) 
    pure_delays=pure_dealy_com(incidence_angle,speed_sound=343,distance_head_center_source=1.5,distance_ear_head=0.09)
    HRTF_fft_fin=HRTF_fft_fin * pure_delays
    HRIR_time_domain = np.fft.ifft(HRTF_fft_fin, axis=-1).real
    full_incidence=np.concatenate((incidence_angle,dist),axis=1)
    for g in range(HRIR_time_domain.shape[0]):
        sota_esti=sf.Sofa('SimpleFreeFieldHRIR')
        sota_esti.Data_IR=HRIR_time_domain[g,:,:,:]
        sota_esti.SourcePosition=full_incidence
        sota_esti.verify()
        sf.write_sofa('/home/PhilipponA/Script_model_HRTF_Sonicom/sofa_hrtf_estimated/'+model_name+'_subject_'+str(g)+'.sofa',sota_esti)

        
    return("uwu")




list_models=os.listdir('/home/PhilipponA/Model_anthro_HRTF/Model_weights/')

device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
for model_name in tqdm(list_models):
        model_path=os.path.join('/home/PhilipponA/Model_anthro_HRTF/Model_weights/',model_name)
        list_models2=os.listdir(model_path)
        list_models2.sort(key=func_sort)
        model=torch.load(os.path.join(model_path,list_models2[-1]),weights_only=False,map_location=device)
        words=model_name.split('_')
        output_type="HRTF_dir"
        val_dataset=HRTF_mesh_dataset(ears_dataset_path='/home/PhilipponA/sonicom_ears_dataset/',HRTF_dataset_path='/databases/sonicom_hrtf_dataset/',type_of_data=words[-1], output_type=output_type,mode=words[0],Train_data=False,Test_data=False,device=device,distorted_photo=True)
        model.eval()
        val_dataloader=torch.utils.data.DataLoader(val_dataset,batch_size=len(val_dataset),shuffle=False)
        incidence_data=val_dataset.incidence_dir
        X_data, Y_data = next(iter(val_dataloader))
        index_inci=list(range(len(incidence_data)))
        HRTF_estimated_all=torch.zeros((X_data.shape[0],793,2,128))
        model_name2=words[-2]+'_'+words[-1]
        if words[-2]!='Manlin-Zhao':
            for i in range(len(incidence_data)):
                inci_gpu=incidence_data.to(device)
                inci_gpu=inci_gpu[i,:]
                inci_gpu=inci_gpu.repeat((X_data.shape[0],1))/360
                out=model(X_data.to(device),inci_gpu.to(device),out_HRTF=True)
                HRTF_estimated_all[:,i,0,:]=out[:,0:128].cpu().detach()
                HRTF_estimated_all[:,i,1,:]=out[:,128:256].cpu().detach()
            HRTF_TO_HRIR_with_phase(HRTF_estimated_all.numpy(),incidence_data.numpy(),model_name2)
        elif words[-2]=='Manlin-Zhao':
            nu=128
            for i in range(128):
                inci_gpu=torch.ones((1),device=device,dtype=torch.long)*i
                inci_gpu=inci_gpu.repeat((X_data.shape[0],1))
                out=model(X_data.to(device),inci_gpu.to(device),out_HRTF=True)
                HRTF_estimated_all[:,:,0,i]=out[:,:,0].cpu().detach()
                HRTF_estimated_all[:,:,1,i]=out[:,:,1].cpu().detach()
            HRTF_TO_HRIR_with_phase(HRTF_estimated_all.numpy(),incidence_data.numpy(),model_name2)


