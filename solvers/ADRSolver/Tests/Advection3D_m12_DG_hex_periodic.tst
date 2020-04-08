<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>3D unsteady DG advection, hexahedra, order 1, P=12,periodic bcs</description>
    <executable>ADRSolver</executable>
    <parameters>Advection3D_m12_DG_hex_periodic.xml</parameters>
    <files>
        <file description="Session File">Advection3D_m12_DG_hex_periodic.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-06">3.56536e-07</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-06">5.38901e-07</value>
        </metric>
    </metrics>
</test>
