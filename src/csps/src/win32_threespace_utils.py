#!/usr/bin/env python2.7

""" This module is a utility module for Windows.
    
    The Win32 ThreeSpace Utils module is a collection of classes, functions,
    structures, and static variables use exclusivly for Windows. All functions
    in this module are used to scan for available ThreeSpace devices on the host
    system and information on them. This module can be used with a system
    running Python 2.5 and newer (including Python 3.x).
"""

__authors__ = [
    '"Chris George" <cgeorge@yeitechnology.com>',
    '"Dan Morrison" <dmorrison@yeitechnology.com>',
]

from threespace_utils import *
import struct
# import serial
# from serial.win32 import ULONG_PTR, is_64bit

import re
import sys
import copy

import ctypes
from ctypes.wintypes import HANDLE
from ctypes.wintypes import BOOL
from ctypes.wintypes import HWND
from ctypes.wintypes import DWORD
from ctypes.wintypes import WORD
from ctypes.wintypes import LONG
from ctypes.wintypes import ULONG
from ctypes.wintypes import LPCSTR
from ctypes.wintypes import HKEY
from ctypes import c_ubyte as BYTE
from ctypes import c_longlong as ULONGLONG
from ctypes import c_wchar as WCHAR
from ctypes import c_ushort as USHORT

from serial.win32 import ULONG_PTR, is_64bit

### Globals ###
NULL = 0
HDEVINFO = ctypes.c_void_p
PCTSTR = ctypes.c_char_p
CHAR = ctypes.c_char
LPDWORD = PDWORD = ctypes.POINTER(DWORD)
LPBYTE = PBYTE = ctypes.c_void_p        # XXX avoids error about types
PHKEY = ctypes.POINTER(HKEY)
ACCESS_MASK = DWORD
REGSAM = ACCESS_MASK
UCHAR = BYTE
BTH_ADDR = ULONGLONG
BLUETOOTH_MAX_NAME_SIZE = 248
HBLUETOOTH_DEVICE_FIND = HANDLE

# Common error enums
ERROR_NO_MORE_ITEMS = 259
ERROR_INVALID_PARAMETER = 87
ERROR_REVISION_MISMATCH = 1306
ERROR_OUTOFMEMORY = 14
ERROR_SUCCESS = 0
ERROR_INVALID_HANDLE = 6
ERROR_MORE_DATA = 234

# hardware enumeration flags
DIGCF_PRESENT = 2
DIGCF_DEVICEINTERFACE = 16
DIGCF_ALLCLASSES = 4
INVALID_HANDLE_VALUE = 0
ERROR_INSUFFICIENT_BUFFER = 122
SPDRP_HARDWAREID = 1
SPDRP_FRIENDLYNAME = 12
ERROR_NO_MORE_ITEMS = 259
DICS_FLAG_GLOBAL = 1
DIREG_DEV = 0x00000001
KEY_READ = 0x20019
REG_SZ = 1
SPDRP_DEVICEDESC = 0
SPDRP_DEVTYPE = 19
SPDRP_DRIVER = 9
SPDRP_ENUMERATOR_NAME  = 0x16
SPDRP_LOCATION_INFORMATION = 0xD
SPDRP_PHYSICAL_DEVICE_OBJECT_NAME = 0xE
SPDRP_MFG = 0xB
SPDRP_SERVICE = 4
SPDRP_CLASS = 7
SPDRP_COMPATIBLEIDS = 2
SPDRP_CLASSGUID = 0x8
SPDRP_ADDRESS = 0x1C

# libraries we use
bthprops = ctypes.windll["bthprops.cpl"]
kernel32 = ctypes.windll["Kernel32.dll"]
setupapi = ctypes.windll.LoadLibrary("setupapi")
advapi32 = ctypes.windll.LoadLibrary("Advapi32")

PortName = b'PortName'

### Classes ###
# COM Port stuctures
class GUID(ctypes.Structure):
    _fields_ = [
        ('Data1', DWORD),
        ('Data2', WORD),
        ('Data3', WORD),
        ('Data4', BYTE * 8),
    ]
    
    def __str__(self):
        return "{%08X-%04X-%04X-%s-%s}" % (
            self.Data1,
            self.Data2,
            self.Data3,
            ''.join(["%02X" % d for d in self.Data4[:2]]),
            ''.join(["%02X" % d for d in self.Data4[2:]]),
        )


class SP_DEVINFO_DATA(ctypes.Structure):
    _fields_ = [
        ('cbSize', DWORD),
        ('ClassGuid', GUID),
        ('DevInst', DWORD),
        ('Reserved', ULONG_PTR),
    ]
    
    def __str__(self):
        return "ClassGuid:%s DevInst:%s" % (self.ClassGuid, self.DevInst)


class SP_DEVICE_INTERFACE_DATA(ctypes.Structure):
    _fields_ = [
        ('cbSize', DWORD),
        ('InterfaceClassGuid', GUID),
        ('Flags', DWORD),
        ('Reserved', ULONG_PTR),
    ]
    
    def __str__(self):
        return "InterfaceClassGuid:%s Flags:%s" % (self.InterfaceClassGuid, self.Flags)


# Bluetooth structures
class BLUETOOTH_DEVICE_SEARCH_PARAMS(ctypes.Structure):
    _fields_ = [
        ('cbSize', DWORD),
        ('fReturnAuthenticated', BOOL),
        ('fReturnRemembered', BOOL),
        ('fReturnUnknown', BOOL),
        ('fReturnConnected', BOOL),
        ('fIssueInquiry', BOOL),
        ('cTimeoutMultiplier', UCHAR),
        ('hRadio', HANDLE),
    ]


class BLUETOOTH_ADDRESS(ctypes.Union):
    _fields_ = [
        ('ullLong', BTH_ADDR),
        ('rgBytes', UCHAR * 6),
    ]
    
    def __str__(self):
        return self.__repr__()
    
    def __repr__(self):
        addr_str = ""
        for i in range(len(self.rgBytes) - 1, -1, -1):
            tmp_str = hex(self.rgBytes[i])[2:].upper()
            if len(tmp_str) < 2:
                tmp_str = "0" + tmp_str
            if i != 0:
                tmp_str += ":"
            addr_str += tmp_str
        return addr_str
    
    def __eq__(self, other):
        if str(self) == str(other):
            return True
        else:
            return False


class SYSTEMTIME(ctypes.Structure):
    _fields_ = [
        ('wYear', WORD),
        ('wMonth', WORD),
        ('wDayOfWeek', WORD),
        ('wDay', WORD),
        ('wHour', WORD),
        ('wMinute', WORD),
        ('wSecond', WORD),
        ('wMilliseconds', WORD),
    ]
    
    def __str__(self):
        month_map = {
            0: "Month_Zero",
            1: "January",
            2: "February",
            3: "March",
            4: "April",
            5: "May",
            6: "June",
            7: "July",
            8: "August",
            9: "September",
            10: "October",
            11: "November",
            12: "December"
        }
        day_of_week_map = {
            0: "Sunday",
            1: "Monday",
            2: "Tuesday",
            3: "Wednesday",
            4: "Thursday",
            5: "Friday",
            6: "Saturday"
        }
        return "%s, %s %d, %d\n%d:%d:%d.%d" % (
            day_of_week_map[self.wDayOfWeek],
            month_map[self.wMonth],
            self.wDay,
            self.wYear,
            self.wHour,
            self.wMinute,
            self.wSecond,
            self.wMilliseconds
        )


class BLUETOOTH_DEVICE_INFO(ctypes.Structure):
    _fields_ = [
        ('cbSize', DWORD),
        ('Address', BLUETOOTH_ADDRESS),
        ('ulClassofDevice', ULONG),
        ('fConnected', BOOL),
        ('fRemembered', BOOL),
        ('fAuthenticated', BOOL),
        ('stLastSeen', SYSTEMTIME),
        ('stLastUsed', SYSTEMTIME),
        ('szName', WCHAR * BLUETOOTH_MAX_NAME_SIZE),
    ]
    
    def __str__(self):
        class_str = hex(self.ulClassofDevice)
        if class_str[-1] == "L":
            class_str = class_str[:-1]
            while len(class_str) < 10:
                class_str = "0x0" + class_str[2:]
        if self.fConnected == 0:
            connected_str = "False"
        else:
            connected_str = "True"
        if self.fRemembered == 0:
            remembered_str = "False"
        else:
            remembered_str = "True"
        if self.fAuthenticated == 0:
            authenticated_str = "False"
        else:
            authenticated_str = "True"
        return (
            "Size: %d\n"            % self.cbSize +
            "Address: %s\n"         % str(self.Address) +
            "Class Of Device: %s\n" % class_str +
            "Connected: %s\n"       % connected_str +
            "Remembered: %s\n"      % remembered_str +
            "Authenticated: %s\n"   % authenticated_str +
            "Last Seen: %s\n"       % str(self.stLastSeen) +
            "Last Used: %s\n"       % str(self.stLastUsed) +
            "Name: %s"              % str(self.szName)
        )

### Helper Functions ###
if sys.version_info >= (3, 0):
    def toLong(number, base=None):
        if base:
            return int(number, base)
        return int(number)
else:
    def toLong(number, base=None):
        if base:
            return long(number, base)
        return long(number)


def _byteBuffer(length):
    return (BYTE * length)()


def _string(buffer):
    s = []
    for c in buffer:
        if c == 0: break
        s.append(chr(c & 0xff)) # "& 0xff": hack to convert signed to unsigned
    return ''.join(s)


def _validHandle(value, func, arguments):
    if value == 0:
        raise ctypes.WinError()
    return value


def _stringToGUID(GUID_string):
    """ Assuming GUID string is formatted as such:
            '{XXXXXXXX-XXXX-XXXX-XXXXXXXXXXXX}'
    """
    return GUID(
        toLong(GUID_string[1:9], 16),
        toLong(GUID_string[10:14], 16),
        toLong(GUID_string[15:19], 16),
        (BYTE * 8)(
            int(GUID_string[20:22], 16),
            int(GUID_string[22:24], 16),
            int(GUID_string[25:27], 16),
            int(GUID_string[27:29], 16),
            int(GUID_string[29:31], 16),
            int(GUID_string[31:33], 16),
            int(GUID_string[33:35], 16),
            int(GUID_string[35:37], 16)
           )
    )


def _stringToBluetoothAddress(address_string):
    """ Assumming address string is formatted as such:
            'XXXXXXXXXXXX'
    """
    tmp_addr = BLUETOOTH_ADDRESS()
    tmp_addr.ullLong = toLong(address_string, 16)
    return tmp_addr


### Structures/Class Pointers ###
PSP_DEVINFO_DATA = ctypes.POINTER(SP_DEVINFO_DATA)
PSP_DEVICE_INTERFACE_DATA = ctypes.POINTER(SP_DEVICE_INTERFACE_DATA)
PSP_DEVICE_INTERFACE_DETAIL_DATA = ctypes.c_void_p

SetupDiDestroyDeviceInfoList = setupapi.SetupDiDestroyDeviceInfoList
SetupDiDestroyDeviceInfoList.argtypes = [HDEVINFO]
SetupDiDestroyDeviceInfoList.restype = BOOL

SetupDiGetClassDevs = setupapi.SetupDiGetClassDevsA
SetupDiGetClassDevs.argtypes = [ctypes.POINTER(GUID), PCTSTR, HWND, DWORD]
SetupDiGetClassDevs.restype = HDEVINFO
SetupDiGetClassDevs.errcheck = _validHandle

SetupDiEnumDeviceInterfaces = setupapi.SetupDiEnumDeviceInterfaces
SetupDiEnumDeviceInterfaces.argtypes = [HDEVINFO, PSP_DEVINFO_DATA, ctypes.POINTER(GUID), DWORD, PSP_DEVICE_INTERFACE_DATA]
SetupDiEnumDeviceInterfaces.restype = BOOL

SetupDiGetDeviceInterfaceDetail = setupapi.SetupDiGetDeviceInterfaceDetailA
SetupDiGetDeviceInterfaceDetail.argtypes = [HDEVINFO, PSP_DEVICE_INTERFACE_DATA, PSP_DEVICE_INTERFACE_DETAIL_DATA, DWORD, PDWORD, PSP_DEVINFO_DATA]
SetupDiGetDeviceInterfaceDetail.restype = BOOL

SetupDiGetDeviceRegistryProperty = setupapi.SetupDiGetDeviceRegistryPropertyA
SetupDiGetDeviceRegistryProperty.argtypes = [HDEVINFO, PSP_DEVINFO_DATA, DWORD, PDWORD, PBYTE, DWORD, PDWORD]
SetupDiGetDeviceRegistryProperty.restype = BOOL

SetupDiOpenDevRegKey = setupapi.SetupDiOpenDevRegKey
SetupDiOpenDevRegKey.argtypes = [HDEVINFO, PSP_DEVINFO_DATA, DWORD, DWORD, DWORD, REGSAM]
SetupDiOpenDevRegKey.restype = HKEY


RegCloseKey = advapi32.RegCloseKey
RegCloseKey.argtypes = [HKEY]
RegCloseKey.restype = LONG

RegQueryValueEx = advapi32.RegQueryValueExA
RegQueryValueEx.argtypes = [HKEY, LPCSTR, LPDWORD, LPDWORD, LPBYTE, LPDWORD]
RegQueryValueEx.restype = LONG

# Used to find 3-Space Sensor devices connected via USB
GUID_DEVINTERFACE_SERENUM_BUS_ENUMERATOR = GUID(toLong(0x4D36E978), 0xE325, 0x11CE, (BYTE * 8)(0xBF, 0xC1, 0x08, 0x00, 0x2B, 0xE1, 0x03, 0x18))

# Used to find Bluetooth and Unknown devices
GUID_DEVINTERFACE_COMPORT = GUID(toLong(0x86E0D1E0), 0x8089, 0x11D0, (BYTE * 8)(0x9C, 0xE4, 0x08, 0x00, 0x3E, 0x30, 0x1F, 0x73))


### Functions ###
def _getBluetoothDevices():
    found_devices = []
    
    ## Create our needed structures
    m_SearchParams = BLUETOOTH_DEVICE_SEARCH_PARAMS()
    m_SearchParams.cbSize = ctypes.sizeof(m_SearchParams)
    m_SearchParams.fReturnAuthenticated = 1 # true
    m_SearchParams.fReturnRemembered = 0 # false
    m_SearchParams.fReturnUnknown = 1 # true
    m_SearchParams.fReturnConnected = 1 # true
    m_SearchParams.fIssueInquiry = 1 # true
    m_SearchParams.cTimeoutMultiplier = 1
    m_SearchParams.hRadio = 0 # Search all available radios
    m_DeviceInfo = BLUETOOTH_DEVICE_INFO()
    m_DeviceInfo.cbSize = ctypes.sizeof(m_DeviceInfo)
    
    device_find_handle = bthprops.BluetoothFindFirstDevice(ctypes.byref(m_SearchParams), ctypes.byref(m_DeviceInfo))
    if device_find_handle == 0:
        # We failed to find a device
        error_code = kernel32.GetLastError()
        # Not sure this is ever returned, but who knows
        if error_code == ERROR_NO_MORE_ITEMS:
            return found_devices
        elif '-d' in sys.argv:
            if error_code == ERROR_INVALID_PARAMETER:
                raise Exception("FindFirstDevice: Either the search params or the device info structure is NULL.")
            elif error_code == ERROR_REVISION_MISMATCH:
                raise Exception("FindFirstDevice: Either the search params or the device info structure is the wrong size.")
            else:
                raise Exception("FindFirstDevice: Unknown function error: %d" % error_code)
    else:
        # We found an initial device
        found_devices.append(copy.deepcopy(m_DeviceInfo))
        while True:
            # Now to find more devices
            found_more_devices = bthprops.BluetoothFindNextDevice(device_find_handle, ctypes.byref(m_DeviceInfo))
            if found_more_devices == 0:
                # We failed to find a device
                error_code = kernel32.GetLastError()
                if error_code == ERROR_NO_MORE_ITEMS:
                    break
                elif '-d' in sys.argv:
                    if error_code == ERROR_INVALID_HANDLE:
                        raise Exception("FindNextDevice: The find handle is NULL.")
                    elif error_code == ERROR_OUTOFMEMORY:
                        raise Exception("FindNextDevice: Out of memory.")
                    else:
                        raise Exception("FindNextDevice: Unknown function error: %d" % error_code)
            else:
                found_devices.append(copy.deepcopy(m_DeviceInfo))
    
    return found_devices


def _yeiGrep(reg_exp):
    for port, desc, hw_id, vid_pid in _yeiComPorts():
        if (re.search(reg_exp, port, re.I) or re.search(reg_exp, desc) or re.search(reg_exp, hw_id)):
            yield port, desc, hw_id, vid_pid


def _yeiComPorts():
    """ This generator scans the device registry for com ports and yields port,
        desc, hw_id
    """
    
    GUID_list = [GUID_DEVINTERFACE_SERENUM_BUS_ENUMERATOR, GUID_DEVINTERFACE_COMPORT]
    ports_yielded = []
    bt_device_list = None
    
    for device_GUID in GUID_list:
        g_hdi = SetupDiGetClassDevs(ctypes.byref(device_GUID), None, NULL, DIGCF_PRESENT|DIGCF_DEVICEINTERFACE)
        for dw_index in range(256):
            friendly_name_string = ""
            did = SP_DEVICE_INTERFACE_DATA()
            did.cbSize = ctypes.sizeof(did)
            
            if not SetupDiEnumDeviceInterfaces(g_hdi, None, ctypes.byref(device_GUID), dw_index, ctypes.byref(did)):
                if ctypes.GetLastError() != ERROR_NO_MORE_ITEMS:
                    if '-d' in sys.argv:
                        raise ctypes.WinError()
                break
            
            dw_needed = DWORD()
            # Get the size
            if not SetupDiGetDeviceInterfaceDetail(g_hdi, ctypes.byref(did), None, 0, ctypes.byref(dw_needed), None):
                # Ignore ERROR_INSUFFICIENT_BUFFER
                if ctypes.GetLastError() != ERROR_INSUFFICIENT_BUFFER:
                    if '-d' in sys.argv:
                        raise ctypes.WinError()
            
            # Allocate buffer
            class SP_DEVICE_INTERFACE_DETAIL_DATA_A(ctypes.Structure):
                _fields_ = [
                    ('cbSize', DWORD),
                    ('DevicePath', CHAR * (dw_needed.value - ctypes.sizeof(DWORD))),
                ]
                
                def __str__(self):
                    return "DevicePath: %s" % self.DevicePath
            
            idd = SP_DEVICE_INTERFACE_DETAIL_DATA_A()
            if is_64bit():
                idd.cbSize = 8
            else:
                idd.cbSize = 5
            
            dev_info = SP_DEVINFO_DATA()
            dev_info.cbSize = ctypes.sizeof(dev_info)
            if not SetupDiGetDeviceInterfaceDetail(g_hdi, ctypes.byref(did), ctypes.byref(idd), dw_needed, None, ctypes.byref(dev_info)):
                if '-d' in sys.argv:
                    raise ctypes.WinError()
            
            # hardware ID
            sz_hardware_id = _byteBuffer(1024)
            if not SetupDiGetDeviceRegistryProperty(g_hdi, ctypes.byref(dev_info), SPDRP_HARDWAREID, None, ctypes.byref(sz_hardware_id), ctypes.sizeof(sz_hardware_id) - 1, None):
                # Ignore ERROR_INSUFFICIENT_BUFFER
                if ctypes.GetLastError() != ERROR_INSUFFICIENT_BUFFER:
                    if '-d' in sys.argv:
                        raise ctypes.WinError()
            
            #Build VID/PID string
            vid_pid_string = ""
            hw_string = _string(sz_hardware_id)
            hw_string = hw_string.upper()
            vid_idx = hw_string.find("VID_")
            pid_idx = hw_string.find("PID_")
            if vid_idx != -1 and pid_idx != -1:
                vid_end = hw_string.find("&", vid_idx + 1)
                vid = hw_string[vid_idx:vid_end]
                pid_end = hw_string.find("&", pid_idx + 1)
                pid = hw_string[pid_idx:pid_end]
                vid_pid_string = vid + "&" + pid
            
            enum_name_buff = _byteBuffer(1024)
            if SetupDiGetDeviceRegistryProperty(g_hdi, ctypes.byref(dev_info), SPDRP_ENUMERATOR_NAME, None, ctypes.byref(enum_name_buff), ctypes.sizeof(enum_name_buff) - 1, None):
                if _string(enum_name_buff).upper() == "BTHENUM":
                    # This is a bluetooth enumerator, we should do further
                    # investigation
                    if bt_device_list is None:
                        bt_device_list = _getBluetoothDevices()
                    
                    device_path_str = idd.DevicePath
                    if type(device_path_str) is bytes:
                        device_path_str = bytes.decode(device_path_str)
                    start_idx = device_path_str.rfind("&") + 1
                    end_idx = start_idx + 12
                    bt_addr_string = device_path_str[start_idx:end_idx]
                    bt_address = _stringToBluetoothAddress(bt_addr_string)
                    if bt_address == _stringToBluetoothAddress("0"):
                        continue
                    connected_dev = None
                    for bt_dev in bt_device_list:
                        if bt_dev.Address == bt_address:
                            connected_dev = bt_dev
                            break
                    if connected_dev is not None:
                        if (str(connected_dev.szName).find("YEI_3SpaceBT") != -1):
                            # The device is a 3-Space Sensor!
                            vid_pid_string = "VID_2476&PID_1060"
                            friendly_name_string = "3 Space Bluetooth over Bluetooth link "
            
            sz_friendly_name = _byteBuffer(1024)
            if not SetupDiGetDeviceRegistryProperty(g_hdi, ctypes.byref(dev_info), SPDRP_FRIENDLYNAME, None, ctypes.byref(sz_friendly_name), ctypes.sizeof(sz_friendly_name) - 1, None):
                # Ignore ERROR_INSUFFICIENT_BUFFER
                if ctypes.GetLastError() != ERROR_INSUFFICIENT_BUFFER:
                    if '-d' in sys.argv:
                        raise IOError("Failed to get details for %s (%s)" % (dev_info, sz_hardware_id.value))
                    port_name = None
            else:
                # The real com port name has to read differently...
                h_key = SetupDiOpenDevRegKey(g_hdi, ctypes.byref(dev_info), DICS_FLAG_GLOBAL, 0, DIREG_DEV, KEY_READ)
                port_name_buffer = _byteBuffer(1024)
                port_name_length = ULONG(ctypes.sizeof(port_name_buffer))
                RegQueryValueEx(h_key, PortName, None, None, ctypes.byref(port_name_buffer), ctypes.byref(port_name_length))
                RegCloseKey(h_key)
                
                # We either use the generated friendly name or our overridden
                # one, with preference to the overridden one.
                if friendly_name_string == "":
                    friendly_name_string = _string(sz_friendly_name)
                else:
                    friendly_name_string += "(" + _string(port_name_buffer) + ")"
                if _string(port_name_buffer) not in ports_yielded:
                    ports_yielded.append(_string(port_name_buffer))
                    yield (_string(port_name_buffer), friendly_name_string, _string(sz_hardware_id), vid_pid_string)
        SetupDiDestroyDeviceInfoList(g_hdi)


def getComPorts(filter=TSS_FIND_ALL):
    """ Queries the system for all available serial COM ports and returns a list
        of them.
        
        Args:
            filter: An interger denoting a flag of what 3-Space Sensors device
                type to be found (default is TSS_FIND_ALL)
        
        Returns:
            A list of all known serial COM ports. Each element of the list is a
                tuple formatted as such:
                    (COM_PORT_NAME, FRIENDLY_NAME, YEI_TECH_DEVICE_TYPE)
            Note:
                YEI_TECH_DEVICE_TYPE will be an empty string if the port's
                    driver's vendor and product IDs do not match any known YEI
                    Techology products.
                Possible YEI_TECH_DEVICE_TYPE strings are:
                    '???' - Unknown
                    'BTL' - Bootloader (No Firmware)
                    'USB' - USB
                    'DNG' - Dongle
                    'WL' - Wireless
                    'EM' - Embedded
                    'DL' - Data-logging
                    'BT' - Bluetooth
    """
    port_list = []
    serial_port_list = _yeiComPorts()
    pid_map = {
        "PID_1000": ("BTL", TSS_FIND_BTL),
        "PID_1010": ("USB", TSS_FIND_USB),
        "PID_1020": ("DNG", TSS_FIND_DNG),
        "PID_1030": ("WL", TSS_FIND_WL),
        "PID_1040": ("EM", TSS_FIND_EM),
        "PID_1050": ("DL", TSS_FIND_DL),
        "PID_1060": ("BT", TSS_FIND_BT)
    }
    for cur_port in serial_port_list:
        hw_string = cur_port[2]
        if cur_port[3] != "":
            vid, pid = cur_port[3].split("&")
            if vid == "VID_2476" and pid in pid_map and pid_map[pid][1] & filter:
                port_list.append(ComInfo(cur_port[0], cur_port[1], pid_map[pid][0]))
                continue
        elif TSS_FIND_UNKNOWN & filter:
            port_list.append(ComInfo(cur_port[0], cur_port[1], "???"))
    return port_list


def _getSoftwareVersionFromPort(serial_port):
    # Figure out whether the current hardware is on "old" or "new" firmware
    serial_port.write(bytearray((0xf7, 0xdf, 0xdf)))
    response = convertString(serial_port.read(9))
    if len(response) == 0:
        # Check and see if in bootloader
        return None
    elif response[:3] == "TSS":
        # Old firmware version remainder
        serial_port.read(9)
        raise Exception("Firmware for device on ( %s ) is out of date for this API. Recommend updating to latest firmware." % serial_port.name)
    
    # Hour-minute remainder
    serial_port.read(3)
    return response


def getDeviceInfoFromComPort(port_name, poll_device=True):
    """ Analyzes a serial COM port of a 3-Space Sensor and returns details about
        the device.
        
        Args:
            port_name: A string representing the name of the serial COM port to
                analyze.
            poll_device: An optional boolean that controls whether the named COM
                port is written to and queried for information about the device.
                If this value is True, please take caution as the COM port's
                device will be written to and may produce undesired effects if
                the device is unknown or not a 3-Space Sensor (default is True)
        
        Returns:
            A list of 5 values describing various details about the COM port's
            device:
                Friendly name,
                3-Space Type,
                3-Space ID,
                3-Space Firmware Version String,
                3-Space Hardware Version String,
                isInBootloader
        
        Raises:
            No explicit exceptions are raised.
    """
    friendly_name = ""
    dev_type = "???"
    dev_serial = 0
    dev_fw_ver = ""
    dev_hw_ver = ""
    in_bootloader = False
    pid_map = {
        "PID_1000": "BTL",
        "PID_1010": "USB",
        "PID_1020": "DNG",
        "PID_1030": "WL",
        "PID_1040": "EM",
        "PID_1050": "DL",
        "PID_1060": "BT"
    }
    matched_ports = _yeiGrep(port_name)

    for cur_port in matched_ports:
        if cur_port[0] == port_name:
            # friendly name
            friendly_name = cur_port[1]
            # device type
            if cur_port[3] != "":
                vid, pid = cur_port[3].split("&")
                if vid == "VID_2476":
                    # The VID matches the YEI vendor ID
                    if pid in pid_map:
                        dev_type = pid_map[pid]
            break
    if poll_device:
        tmp_port = None
        try:
            tmp_port = serial.Serial(port_name, timeout=0.1, baudrate=115200)
            
            if dev_type == "BT":
                tmp_port.timeout = 5.0
        except:
            tmp_port = None
        if tmp_port is not None:
            # Try to get the serial, if it fails try to see if in bootloader
            tmp_port.write(bytearray((0xf7, 0xed, 0xed)))
            response = tmp_port.read(4)
            if len(response) == 4:
                dev_serial = "{0:08X}".format(struct.unpack('>I', response)[0])
                # Get the version strings (and device type if the
                # previous method did not resolve it)
                software_version = _getSoftwareVersionFromPort(tmp_port)
                if software_version is not None:
                    # This is in fact a 3-Space sensor
                    dev_fw_ver = software_version
                    tmp_port.write(bytearray((0xf7, 0xe6, 0xe6)))
                    hardware_version = convertString(tmp_port.read(32))
                    dev_hw_ver = hardware_version
                    if dev_type == "???":
                        dev_type = hardware_version[4:-8].strip()
                else:
                    tmp_port.write(bytearray((0x3f,))) # this is ascii '?'
                    response = convertString(tmp_port.read(2))
                    if response:
                        if response == "OK":
                            in_bootloader = True
                            dev_type = "BTL"
                    else:
                        raise Exception("Either device on( %s ) is not a 3-Space Sensor or the firmware is out of date for this API and recommend updating to latest firmware." % port_name)
            
            tmp_port.close()
    return SensorInfo(
        friendly_name,
        dev_type,
        dev_serial,
        dev_fw_ver,
        dev_hw_ver,
        in_bootloader
    )