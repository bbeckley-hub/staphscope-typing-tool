#!/bin/bash

if [ "$1" == "path" ]
then
    echo "Testing stxtyper command in your \$PATH"
    which stxtyper
    STXTYPER="stxtyper --debug"
else
    echo "Testing ./stxtyper"
    echo "  To test stxtyper in your path run 'test_stxtyper.sh path'"
    STXTYPER="./stxtyper  --debug"
fi

if [ ! -e "test/basic.expected" ]
then
    echo "test/basic.expected not found, downloading new test data"
    echo "from https://raw.githubusercontent.com/ncbi/stxtyper/main/"
    mkdir test
    pushd test
        curl --location --silent -O https://raw.githubusercontent.com/ncbi/stxtyper/main/test/amrfinder_integration.expected \
             -O https://raw.githubusercontent.com/ncbi/stxtyper/main/test/amrfinder_integration.fa \
             -O https://raw.githubusercontent.com/ncbi/stxtyper/main/test/amrfinder_integration2.expected \
             -O https://raw.githubusercontent.com/ncbi/stxtyper/main/test/amrfinder_integration2.fa \
             -O https://raw.githubusercontent.com/ncbi/stxtyper/main/test/basic.expected \
             -O https://raw.githubusercontent.com/ncbi/stxtyper/main/test/basic.fa \
             -O https://raw.githubusercontent.com/ncbi/stxtyper/main/test/basic.nuc_out.expected \
             -O https://raw.githubusercontent.com/ncbi/stxtyper/main/test/cases.expected \
             -O https://raw.githubusercontent.com/ncbi/stxtyper/main/test/cases.fa \
             -O https://raw.githubusercontent.com/ncbi/stxtyper/main/test/synthetics.expected \
             -O https://raw.githubusercontent.com/ncbi/stxtyper/main/test/synthetics.fa \
             -O https://raw.githubusercontent.com/ncbi/stxtyper/main/test/virulence_ecoli.expected \
             -O https://raw.githubusercontent.com/ncbi/stxtyper/main/test/virulence_ecoli.fa
    popd
fi

# globals updated by function test_input_file
FAILURES=0
TESTS=0
TEST_TEXT=''

# echo "TERM=$TERM"

# some color macros
if [ "$TERM" == "" ] || [ "$TERM" == "dumb" ] || [ ! -t 1 ]
then
    green='' # no colors
    red=''
    bold=''
    reset=''
else
    green=`tput setaf 2`  # Set green foreground color (code 2)
    red=`tput setaf 1`    # Set red foreground color (code 1)
    bold=`tput bold`      # Set bold
    reset=`tput sgr0`     # Reset color to default
fi

function test_input_file {
    local test_base="$1"
    local options="$2"

    TESTS=$(( $TESTS + 1 ))

    if ! $STXTYPER $options -n "test/$test_base.fa" > "test/$test_base.got"
    then
        echo "${red}not ok: $STXTYPER returned a non-zero exit value indicating a failure of the software${reset}"
        echo "#  $STXTYPER $options -n test/$test_base.fa > test/$test_base.got"
        TEST_TEXT="$TEST_TEXT"$'\n'"${red}Failed $test_base${reset}"
        return 1
    else
        if ! diff -q "test/$test_base.expected" "test/$test_base.got"
        then
            echo "${red}not ok: $STXTYPER returned output different from expected${reset}"
            echo "#  $STXTYPER $options -n test/$test_base.fa > test/$test_base.got"
            echo "# diff test/$test_base.expected test/$test_base.got"
            diff "test/$test_base.expected" "test/$test_base.got"
            echo "#  To approve run:"
            echo "#     mv test/$test_base.got test/$test_base.expected "
            TEST_TEXT="$TEST_TEXT"$'\n'"${red}Failed $test_base${reset}"
            return 1
        else
            echo "${green}ok:${reset} test/$test_base.fa"
            return 0
        fi
    fi
}


test_input_file 'basic' '--nucleotide_output test/basic.nuc_out.got'
FAILURES=$(( $? + $FAILURES ))
# test --nucleotide_output option
TESTS=$(( $TESTS + 1 ))
if ! diff -q "test/basic.nuc_out.expected" "test/basic.nuc_out.got"
then
    echo "${red}not ok: $STXTYPER returned --nucleotide_output file different from expected${reset}"
    echo "#  $STXTYPER --nucleotide_output test/basic.nuc_out.got -n test/basic.fa > test/basic.got"
    echo "# diff test/basic.nuc_out.expected test/basic.nuc_out.got"
    diff "test/basic.nuc_out.expected" "test/basic.nuc_out.got"
    echo "#  To approve run:"
    echo "#     mv test/basic.nuc_out.got test/basic.nuc_out.expected "
    TEST_TEXT="$TEST_TEXT"$'\n'"${red}Failed basic --nucleotide_output test${reset}"
    FAILURES=$(( $FAILURES + 1 ))
else 
    echo "${green}ok:${reset} --nucleotide_output test/basic.nuc_out.got options worked"
fi

test_input_file 'synthetics'
FAILURES=$(( $? + $FAILURES ))

test_input_file 'virulence_ecoli' '--threads 4'
FAILURES=$(( $? + $FAILURES ))

test_input_file 'cases'
FAILURES=$(( $? + $FAILURES ))

test_input_file 'amrfinder_integration' '--amrfinder'
FAILURES=$(( $? + $FAILURES ))

test_input_file 'amrfinder_integration2' '--amrfinder --print_node' 
FAILURES=$(( $? + $FAILURES ))

echo "Done."
echo "$TEST_TEXT"
echo ""
if [ "$FAILURES" -gt 0 ]
then
    PASSED=$(( $TESTS - $FAILURES ))
    echo "${red}not ok overall: $FAILURES out of $TESTS stxtyper tests failed${reset}"
    exit 1
else
    echo "${green}ok: all $TESTS stxtyper tests passed${reset}"
    echo "${green}${bold}Success!${reset}"
fi
