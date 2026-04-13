import tempfile
import shutil
import os
import glob


def move_all_contents(src_dir, dst_dir):
    # Ensure the source directory exists
    if not os.path.exists(src_dir):
        raise FileNotFoundError(f"Source directory '{src_dir}' does not exist.")

    # Ensure the destination directory exists
    if not os.path.exists(dst_dir):
        os.makedirs(dst_dir)

    # Iterate over all items in the source directory
    for item in os.listdir(src_dir):
        src_item = os.path.join(src_dir, item)
        dst_item = os.path.join(dst_dir, item)

        # Move the item to the destination directory
        shutil.move(src_item, dst_item)


class TemporaryFolderOverride:
    def __init__(self, target_folder):
        self.target_folder_exists = os.path.exists(target_folder)
        self.target_folder = target_folder
        self.backup_folder = None

    def __enter__(self):
        # Create a temporary directory
        if self.target_folder_exists:
            self.backup_folder = tempfile.mkdtemp()
            move_all_contents(self.target_folder, self.backup_folder)
        else:
            # Create a new empty folder in place of the original
            os.mkdir(self.target_folder)

    def __exit__(self, exc_type, exc_value, traceback):
        try:
            # Empty the target folder
            for path in glob.glob(os.path.join(self.target_folder, '*')):
                shutil.rmtree(path) if os.path.isdir(path) and not os.path.islink(path) else os.unlink(path)
            if self.target_folder_exists:
                move_all_contents(self.backup_folder, self.target_folder)
            # If the target folder is not a mount point and it did not exist before
            elif self.target_folder_exists and not os.path.ismount(self.target_folder):
                shutil.rmtree(self.target_folder)

        finally:
            if self.backup_folder is not None:
                shutil.rmtree(self.backup_folder, ignore_errors=True)
