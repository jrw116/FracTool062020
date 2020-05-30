from abaqusGui import getAFXApp, Activator, AFXMode, afxCreatePNGIcon
from abaqusConstants import ALL
import os
thisPath = os.path.abspath(__file__)
thisDir = os.path.dirname(thisPath)

toolset = getAFXApp().getAFXMainWindow().getPluginToolset()
toolset.registerGuiMenuButton(
    buttonText='FracTool',
    object=Activator(os.path.join(thisDir, 'plugin3DB.py')),
    kernelInitString='import EDC2',
    messageId=AFXMode.ID_ACTIVATE,
    icon=afxCreatePNGIcon(os.path.join(thisDir,'logo.png')),
    applicableModules=ALL,
    version='N/A',
    author='N/A',
    description='N/A',
    helpUrl='N/A'
)
