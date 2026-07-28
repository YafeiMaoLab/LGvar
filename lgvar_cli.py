import os
import subprocess
import sys


def main():
    package_dir = os.path.dirname(os.path.abspath(__file__))
    lgvar_script = os.path.join(package_dir, "LGVAR")

    if os.path.exists(lgvar_script):
        os.chmod(lgvar_script, 0o755)
        cmd = [lgvar_script] + sys.argv[1:]
        try:
            result = subprocess.run(cmd, cwd=package_dir)
            sys.exit(result.returncode)
        except Exception as e:
            print(f"[Error] Failed to execute LGVAR: {e}")
            sys.exit(1)
    else:
        print(f"[Error] Could not locate 'LGVAR' in {package_dir}")
        sys.exit(1)


if __name__ == "__main__":
    main()
