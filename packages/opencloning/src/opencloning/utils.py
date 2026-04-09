import tempfile
import shutil
import os


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
        # Remove the temporary folder
        shutil.rmtree(self.target_folder)

        if self.target_folder_exists:
            os.mkdir(self.target_folder)
            move_all_contents(self.backup_folder, self.target_folder)
