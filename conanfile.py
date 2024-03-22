from conan import ConanFile
from conan.tools.cmake import CMake
from conan.tools.files import copy
import os

class Project(ConanFile):

    # all project are the same:
    python_requires = "project_base/1.0"
    python_requires_extend = "project_base.ProjectBase"

    def init(self):
        base = self.python_requires["project_base"].module.ProjectBase
        self.settings = base.settings
        self.options.update(base.options, base.default_options)
        self.revision_mode = base.revision_mode

    def set_name(self):
        self.name = "opengv"

    # difference between project:
    def requirements(self):
        self.requires("eigen/3.3.7", transitive_libs=True)
        self.requires("glog/0.6.0", transitive_libs=True)
        
    def package_info(self):
        self.cpp_info.set_property("cmake_file_name", "opengv")
        self.cpp_info.set_property("cmake_target_name", "opengv::opengv")
        self.cpp_info.includedirs = ["include"]
        self.cpp_info.libs = ["camera_models"]
