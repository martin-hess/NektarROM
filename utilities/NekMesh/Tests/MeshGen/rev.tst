<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description>High AoA naca0012 wing with bl</description>
    <executable>NekMesh</executable>
    <parameters>-m jac:list rev.mcf rev.xml:xml:test</parameters>
    <files>
        <file description="Input File">rev.mcf</file>
        <file description="Input File 2">rev-rotated.STEP</file>
    </files>
    <metrics>
        <metric type="regex" id="1">
            <regex>^Total negative Jacobians: (\d+)</regex>
            <matches>
                <match>
                    <field id="0">0</field>
                </match>
            </matches>
        </metric>
    </metrics>
</test>
