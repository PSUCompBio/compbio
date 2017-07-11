install(FILES csps.py threespace_api.py threespace_utils.py win32_threespace_utils.py DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

install (CODE "execute_process (COMMAND pyinstaller --noconfirm \
							--clean \
							--log-level=ERROR \
							--noupx \
							--onedir \
							--icon="./CSPS_HR.ico" \
							--add-data=./BS:./BS \
							--add-data=./Simulation:./Simulation \
							--add-data=./CSPS_HR.ico:./ \
							csps.py WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})")
