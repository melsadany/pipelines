import whisperx
import csv
import argparse
from pathlib import Path
import pandas as pd
from multiprocessing import Pool, cpu_count
import os
import torch  # Explicitly import torch


def process_single_file(args):
    """Process a single audio file"""
    audio_path, output_path, model_name, device_name, compute_type = args
    # Skip if output already exists
    if Path(output_path).exists():
        print(f"Skipping {audio_path}, output already exists")
        return
    try:
        # Load model inside the worker process (with proper torch context)
        import torch
        device = device_name
        model = whisperx.load_model(model_name, device, compute_type=compute_type)
        # Load and transcribe audio
        audio = whisperx.load_audio(audio_path)
        result = model.transcribe(audio, language="en")
        # Load alignment model and align
        align_model, metadata = whisperx.load_align_model(language_code="en", device=device)
        result_aligned = whisperx.align(result["segments"], align_model, metadata, audio, device)
        # Save as TSV
        with open(output_path, "w", newline='') as f:
            writer = csv.writer(f, delimiter='\t')
            writer.writerow(["start", "end", "word"])
            for segment in result_aligned["segments"]:
                for word in segment["words"]:
                    writer.writerow([
                        word.get("start", ""),
                        word.get("end", ""),
                        word.get("word", "")
                    ])
        print(f"Processed: {audio_path}")
    except Exception as e:
        print(f"Error processing {audio_path}: {str(e)}")


def main():
    parser = argparse.ArgumentParser(description="Transcribe and align multiple audio files with WhisperX.")
    parser.add_argument("csv_file", type=str, help="Path to CSV file with columns: audio_path, output_path")
    parser.add_argument("--workers", type=int, default=None, help="Number of parallel workers (default: CPU count)")
    parser.add_argument("--model", type=str, default="large-v3", help="Whisper model to use")
    parser.add_argument("--compute-type", type=str, default="float32", help="Compute type")
    args = parser.parse_args()

    # Load CSV
    df = pd.read_csv(args.csv_file)
    
    # Determine device
    device = "cuda" if torch.cuda.is_available() else "cpu"
    print(f"Using device: {device}")
    
    # Prepare arguments for each file (pass model name instead of model object)
    tasks = [(row.audio_path, row.output_path, args.model, device, args.compute_type) 
             for _, row in df.iterrows()]
    
    # Process in parallel
    num_workers = args.workers or min(cpu_count(), len(tasks))
    print(f"Processing {len(tasks)} files with {num_workers} workers...")
    
    with Pool(num_workers) as pool:
        pool.map(process_single_file, tasks)
    
    print("All files processed!")

if __name__ == "__main__":
    import torch  # Import torch at module level
    main()
