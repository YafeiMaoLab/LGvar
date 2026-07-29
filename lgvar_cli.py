import os
import subprocess
import sys
from datetime import datetime

def log_info(message):
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S,%f")[:-3]
    print(f"{timestamp} - INFO - {message}")

def log_error(message):
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S,%f")[:-3]
    print(f"{timestamp} - ERROR - {message}", file=sys.stderr)

def main():
    package_dir = os.path.dirname(os.path.abspath(__file__))
    
    if len(sys.argv) > 1 and sys.argv[1] == "setup":
        setup_script = os.path.join(package_dir, "setup.sh")
        if os.path.exists(setup_script):
            os.chmod(setup_script, 0o755)
            
            if "--spawn" in sys.argv or "-s" in sys.argv:
                log_info("Spawning a new bash subshell with LGVAR environment loaded...")
                bash_cmd = (
                    f"test -f ~/.bashrc && source ~/.bashrc; "
                    f"source '{setup_script}' && exec bash -i"
                )
                
                subprocess.run(["/bin/bash", "-c", bash_cmd])
                sys.exit(0)
            else:
                log_info("LGVAR Environment Setup")
                log_info("To load LGVAR environment variables into your current shell, run:")
                log_info(f"  source {setup_script}")
                log_info("Or spawn a new configured shell directly using:")
                log_info("  LGVAR setup --spawn")
                sys.exit(0)
        else:
            log_error(f"Could not locate 'setup.sh' in {package_dir}")
            sys.exit(1)

    lgvar_script = os.path.join(package_dir, "LGVAR")
    if os.path.exists(lgvar_script):
        os.chmod(lgvar_script, 0o755)
        cmd = [lgvar_script] + sys.argv[1:]
        try:
            result = subprocess.run(cmd, cwd=package_dir)
            sys.exit(result.returncode)
        except Exception as e:
            log_error(f"Failed to execute LGVAR: {e}")
            sys.exit(1)
    else:
        log_error(f"Could not locate 'LGVAR' in {package_dir}")
        sys.exit(1)

if __name__ == "__main__":
    main()
